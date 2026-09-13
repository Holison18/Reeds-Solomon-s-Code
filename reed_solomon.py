"""
This code is an implementation of Reed-Solomon error correction over GF(256), suitable for QR codes and similar applications. 
It provides functions for encoding messages with parity bytes, as well as decoding and correcting errors and erasures in received messages.
"""

import random
from typing import List, Optional

#GF(256) arithmetic

FIELD_SIZE = 256
PRIM_POLY = 0x11D  # x^8 + x^4 + x^3 + x^2 + 1 -- the field polynomial QR codes use


def _build_gf_tables():
    """Build log/antilog tables for GF(256) using generator alpha = 2."""
    exp = [0] * (2 * FIELD_SIZE)   # oversized so mul/div never need a modulo
    log = [0] * FIELD_SIZE
    x = 1
    for i in range(FIELD_SIZE - 1):
        exp[i] = x
        log[x] = i
        x <<= 1
        if x & FIELD_SIZE:      # overflowed 8 bits -> reduce mod PRIM_POLY
            x ^= PRIM_POLY
    for i in range(FIELD_SIZE - 1, len(exp)):
        exp[i] = exp[i - (FIELD_SIZE - 1)]
    return exp, log


GF_EXP, GF_LOG = _build_gf_tables()


def gf_mul(a: int, b: int) -> int:
    if a == 0 or b == 0:
        return 0
    return GF_EXP[GF_LOG[a] + GF_LOG[b]]


def gf_div(a: int, b: int) -> int:
    if b == 0:
        raise ZeroDivisionError("division by zero in GF(256)")
    if a == 0:
        return 0
    return GF_EXP[(GF_LOG[a] - GF_LOG[b]) % (FIELD_SIZE - 1)]


def gf_pow(a: int, power: int) -> int:
    return GF_EXP[(GF_LOG[a] * power) % (FIELD_SIZE - 1)]


def gf_inverse(a: int) -> int:
    return GF_EXP[(FIELD_SIZE - 1) - GF_LOG[a]]


# Polynomials over GF(256)


def gf_poly_scale(p: List[int], scalar: int) -> List[int]:
    return [gf_mul(c, scalar) for c in p]


def gf_poly_add(p: List[int], q: List[int]) -> List[int]:
    r = [0] * max(len(p), len(q))
    for i, c in enumerate(p):
        r[i + len(r) - len(p)] = c
    for i, c in enumerate(q):
        r[i + len(r) - len(q)] ^= c
    return r


def gf_poly_mul(p: List[int], q: List[int]) -> List[int]:
    r = [0] * (len(p) + len(q) - 1)
    for j, qc in enumerate(q):
        if qc == 0:
            continue
        for i, pc in enumerate(p):
            r[i + j] ^= gf_mul(pc, qc)
    return r


def gf_poly_eval(p: List[int], x: int) -> int:
    """Evaluate p(x) using Horner's method."""
    y = p[0]
    for c in p[1:]:
        y = gf_mul(y, x) ^ c
    return y


# Encoding


def rs_generator_poly(nsym: int) -> List[int]:
    """g(x) = (x - a^0)(x - a^1)...(x - a^(nsym-1))"""
    g = [1]
    for i in range(nsym):
        g = gf_poly_mul(g, [1, GF_EXP[i]])
    return g


def rs_encode_msg(msg_in: List[int], nsym: int) -> List[int]:
    """Systematic RS encode: returns msg_in followed by nsym parity bytes."""
    gen = rs_generator_poly(nsym)
    msg_out = list(msg_in) + [0] * nsym
    for i in range(len(msg_in)):
        coef = msg_out[i]
        if coef != 0:
            for j in range(len(gen)):
                msg_out[i + j] ^= gf_mul(gen[j], coef)
    msg_out[:len(msg_in)] = msg_in  # restore data bytes (long-division trashes them)
    return msg_out


# Decoding

def rs_calc_syndromes(msg: List[int], nsym: int) -> List[int]:
    """Zero everywhere iff msg is a valid codeword."""
    return [gf_poly_eval(msg, GF_EXP[i]) for i in range(nsym)]


def rs_find_error_locator(synd: List[int], nsym: int,
                           erase_loc: Optional[List[int]] = None,
                           erase_count: int = 0) -> List[int]:
    """Berlekamp-Massey algorithm. If erasures are known, seed with their
    locator polynomial so the algorithm only has to find the remaining
    (unlocated) errors."""
    if erase_loc is not None:
        err_loc = list(erase_loc)
        old_loc = list(erase_loc)
    else:
        err_loc = [1]
        old_loc = [1]

    for i in range(nsym - erase_count):
        k = (erase_count + i) if erase_loc is not None else i
        delta = synd[k]
        for j in range(1, len(err_loc)):
            delta ^= gf_mul(err_loc[-(j + 1)], synd[k - j])
        old_loc = old_loc + [0]
        if delta != 0:
            if len(old_loc) > len(err_loc):
                new_loc = gf_poly_scale(old_loc, delta)
                old_loc = gf_poly_scale(err_loc, gf_inverse(delta))
                err_loc = new_loc
            err_loc = gf_poly_add(err_loc, gf_poly_scale(old_loc, delta))

    while err_loc and err_loc[0] == 0:
        err_loc.pop(0)
    errs = len(err_loc) - 1
    if (errs - erase_count) * 2 + erase_count > nsym:
        raise ValueError("too many errors/erasures to correct")
    return err_loc


def rs_find_errors(err_loc: List[int], nmess: int) -> List[int]:
    """Chien search: trial-substitute every position to find the roots of
    the error locator polynomial, i.e. the error positions."""
    errs = len(err_loc) - 1
    err_pos = []
    for i in range(nmess):
        if gf_poly_eval(err_loc, GF_EXP[255 - i]) == 0:
            err_pos.append(nmess - 1 - i)
    if len(err_pos) != errs:
        raise ValueError("error locator has wrong number of roots -- too many errors")
    return err_pos


def _poly_mul_ascending(p: List[int], q: List[int]) -> List[int]:
    """Multiply two polynomials given in ASCENDING order (index i = coefficient
    of x^i). Only used internally by the Forney step below, which is easier
    to state correctly in this convention than in the descending one used
    everywhere else in this module."""
    r = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a == 0:
            continue
        for j, b in enumerate(q):
            r[i + j] ^= gf_mul(a, b)
    return r


def _poly_eval_ascending(p: List[int], x: int) -> int:
    """Evaluate an ascending-order polynomial at x."""
    result = 0
    x_pow = 1
    for c in p:
        result ^= gf_mul(c, x_pow)
        x_pow = gf_mul(x_pow, x)
    return result


def rs_correct_errata(msg: List[int], synd: List[int], pos: List[int]) -> None:
    """Forney algorithm: given known error/erasure positions, compute the
    correction value for each and apply it to msg in place.

    Uses the standard key-equation approach: build the error locator
    Lambda(x) = prod(1 - X_k x) over the error locations X_k, derive the
    error evaluator Omega(x) = [S(x) Lambda(x)] mod x^nsym, then recover
    each error magnitude as Omega(X_k^-1) / prod_{j!=k}(1 - X_j X_k^-1).
    """
    n = len(msg)
    locations = [GF_EXP[n - 1 - p] for p in pos]  # X_k = alpha^(degree of position p)

    lam = [1]  # Lambda(x), ascending order
    for x_k in locations:
        new_lam = [0] * (len(lam) + 1)
        for i, c in enumerate(lam):
            new_lam[i] ^= c
            new_lam[i + 1] ^= gf_mul(c, x_k)
        lam = new_lam

    nsym = len(synd)
    omega = _poly_mul_ascending(synd, lam)[:nsym]  # Omega(x) mod x^nsym

    for k, x_k in enumerate(locations):
        x_k_inv = gf_inverse(x_k)
        numerator = _poly_eval_ascending(omega, x_k_inv)
        denominator = 1
        for j, x_j in enumerate(locations):
            if j != k:
                denominator = gf_mul(denominator, 1 ^ gf_mul(x_k_inv, x_j))
        msg[pos[k]] ^= gf_div(numerator, denominator)


def rs_correct_msg(msg_in: List[int], nsym: int,
                    erase_pos: Optional[List[int]] = None) -> List[int]:
    """Full error-and-erasure correcting decode. Returns the corrected
    message (data + parity), or raises ValueError if it can't be fixed."""
    if len(msg_in) > 255:
        raise ValueError("message too long for GF(256)-based RS")
    msg_out = list(msg_in)
    erase_pos = list(erase_pos) if erase_pos else []
    if len(erase_pos) > nsym:
        raise ValueError("too many erasures to correct")

    synd = rs_calc_syndromes(msg_out, nsym)
    if max(synd) == 0:
        return msg_out  # no errors at all

    erase_loc = None
    if erase_pos:
        erase_loc = [1]
        for p in erase_pos:
            x = GF_EXP[len(msg_out) - 1 - p]
            erase_loc = gf_poly_mul(erase_loc, [x, 1])

    err_loc = rs_find_error_locator(synd, nsym, erase_loc=erase_loc, erase_count=len(erase_pos))
    # err_loc's roots already include the known erasure positions (it was
    # seeded with erase_loc above), so the Chien search below finds the
    # complete set of errata positions -- erasures and newly-found errors
    # together -- with nothing left to concatenate.
    errata_pos = rs_find_errors(err_loc, len(msg_out))

    rs_correct_errata(msg_out, synd, errata_pos)
    synd = rs_calc_syndromes(msg_out, nsym)
    if max(synd) != 0:
        raise ValueError("decoding failed -- message still invalid after correction")
    return msg_out


# Demonstration with text messages (e.g., "He is a boy")


def demo_text_message(text: str = "He is a boy", nsym: int = 8):
    """Demonstrates how Reed-Solomon encodes a plain text string, withstands
    corruptions (errors and erasures), and prevents information loss."""
    print("=" * 65)
    print(f"Reed-Solomon Text Demonstration: \"{text}\"")
    print("=" * 65)

    data_bytes = list(text.encode("utf-8"))
    max_errors = nsym // 2
    max_erasures = nsym

    print(f"Original Text:      '{text}'")
    print(f"Data Bytes:         {data_bytes} (length: {len(data_bytes)})")
    print(f"Parity Symbols:     {nsym} symbols")
    print(f"Protection:         Can correct up to {max_errors} unknown error(s) OR {max_erasures} known erasure(s)\n")

    # 1. Encode message
    codeword = rs_encode_msg(data_bytes, nsym)
    parity_bytes = codeword[len(data_bytes):]
    print(f"Encoded Codeword:   {codeword}")
    print(f"Parity Bytes:       {parity_bytes}\n")

    # 2. Unknown Errors Simulation (noise / data corruption)
    damaged = list(codeword)
    # Corrupt characters in the text (e.g. index 1: 'e' -> '?', index 4: 's' -> 'X', index 9: 'o' -> '0')
    error_positions = [1, 4, 9][:max_errors]
    corrupt_chars = ['?', 'X', '0']
    for idx, pos in enumerate(error_positions):
        damaged[pos] = ord(corrupt_chars[idx])

    corrupted_text = bytes(damaged[:len(data_bytes)]).decode("utf-8", errors="replace")
    print("--- Simulation 1: Unknown Errors (Corruption) ---")
    print(f"Corrupted Text:     '{corrupted_text}' (errors introduced at indices {error_positions})")
    print(f"Corrupted Bytes:    {damaged}")

    # Correct without knowing error positions
    recovered_codeword = rs_correct_msg(damaged, nsym)
    recovered_data = recovered_codeword[:len(data_bytes)]
    recovered_text = bytes(recovered_data).decode("utf-8", errors="replace")

    print(f"Recovered Text:     '{recovered_text}'")
    assert recovered_text == text, "Failed to recover from errors!"
    print("Result:             SUCCESS - Original message restored with zero information loss!\n")

    # 3. Known Erasures Simulation (packet loss / missing bytes)
    erased = list(codeword)
    erasure_positions = [0, 3, 6, 8]  # positions where characters were completely lost
    for p in erasure_positions:
        erased[p] = 0

    erased_display = "".join(chr(erased[i]) if i not in erasure_positions else "_" for i in range(len(data_bytes)))
    print("--- Simulation 2: Erasures (Lost / Missing Bytes) ---")
    print(f"Received Text:      '{erased_display}' (bytes lost at indices {erasure_positions})")

    recovered_erased_codeword = rs_correct_msg(erased, nsym, erase_pos=erasure_positions)
    recovered_erased_data = recovered_erased_codeword[:len(data_bytes)]
    recovered_erased_text = bytes(recovered_erased_data).decode("utf-8", errors="replace")

    print(f"Recovered Text:     '{recovered_erased_text}'")
    assert recovered_erased_text == text, "Failed to recover from erasures!"
    print("Result:             SUCCESS - Missing characters completely reconstructed from parity!\n")
    print("=" * 65 + "\n")


# Testing
# using the QR sample message and its known parity bytes.


if __name__ == "__main__":
    data = [0x40, 0xd2, 0x75, 0x47, 0x76, 0x17, 0x32, 0x06,
            0x27, 0x26, 0x96, 0xc6, 0xc6, 0x96, 0x70, 0xec]
    expected_parity = [0xbc, 0x2a, 0x90, 0x13, 0x6b, 0xaf, 0xef, 0xfd, 0x4b, 0xe0]
    nsym = 10

    encoded = rs_encode_msg(data, nsym)
    computed_parity = encoded[len(data):]
    print("computed parity:", [hex(b) for b in computed_parity])
    assert computed_parity == expected_parity, "encoder does not match the article's worked example"
    print("encode: OK, matches the article's worked example\n")

    # syndromes of a clean codeword must be all zero
    assert max(rs_calc_syndromes(encoded, nsym)) == 0
    print("syndrome check on clean codeword: OK (all zero)\n")

    # damage up to nsym//2 = 5 bytes at unknown positions and recover them
    random.seed(42)
    damaged = list(encoded)
    max_errors = nsym // 2
    error_positions = random.sample(range(len(damaged)), max_errors)
    for p in error_positions:
        damaged[p] ^= random.randint(1, 255)
    print(f"introduced {max_errors} errors at positions {sorted(error_positions)}")

    recovered = rs_correct_msg(damaged, nsym)
    assert recovered == encoded, "failed to recover the original codeword"
    print("error correction (unknown positions): OK, message fully recovered\n")

    # erasure correction: up to nsym erasures if positions are known
    erased = list(encoded)
    erase_positions = random.sample(range(len(erased)), nsym)
    for p in erase_positions:
        erased[p] = -1
    erased_clean = [0 if v < 0 else v for v in erased]
    recovered2 = rs_correct_msg(erased_clean, nsym, erase_pos=erase_positions)
    assert recovered2 == encoded, "failed to recover from erasures"
    print(f"erasure correction ({nsym} known positions): OK, message fully recovered\n")

    # Run the text demonstration
    demo_text_message("He is a boy", nsym=8)

