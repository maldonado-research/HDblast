"""Independent MPF Fourier accumulation with fixed 64-sample refresh blocks."""
from mpmath.libmp import mpf_add, mpf_sub, mpf_mul, mpf_cos, mpf_sin, mpf_neg, fzero, round_nearest

REFRESH_SAMPLES = 64


def accumulate_node(z_initial, omega, k_squared, first_offset, time_step, sample_count,
                    sums, precision, factor=REFRESH_SAMPLES):
    """Accumulate Re(Z), Im(Omega Z), Re(k² Z) without a per-sample transcendental."""
    add, mul, sub = mpf_add, mpf_mul, mpf_sub
    rnd = round_nearest
    step_angle = mul(omega, time_step, precision, rnd)
    cr = mpf_cos(step_angle, precision, rnd)
    ci = mpf_sin(step_angle, precision, rnd)
    offset = first_offset
    for block_start in range(0, sample_count, factor):
        angle = mul(omega, offset, precision, rnd)
        er = mpf_cos(angle, precision, rnd)
        ei = mpf_sin(angle, precision, rnd)
        zr = sub(mul(z_initial[0], er, precision, rnd), mul(z_initial[1], ei, precision, rnd), precision, rnd)
        zi = add(mul(z_initial[0], ei, precision, rnd), mul(z_initial[1], er, precision, rnd), precision, rnd)
        block_end = min(sample_count, block_start + factor)
        for index in range(block_start, block_end):
            sums[0][index] = add(sums[0][index], zr, precision, rnd)
            sums[1][index] = add(sums[1][index], mul(omega, zi, precision, rnd), precision, rnd)
            sums[2][index] = add(sums[2][index], mul(k_squared, zr, precision, rnd), precision, rnd)
            nr = sub(mul(zr, cr, precision, rnd), mul(zi, ci, precision, rnd), precision, rnd)
            zi = add(mul(zr, ci, precision, rnd), mul(zi, cr, precision, rnd), precision, rnd)
            zr = nr
        for _ in range(block_end - block_start):
            offset = add(offset, time_step, precision, rnd)
    return sums


def empty_sums(samples):
    return [[fzero for _ in range(samples)] for _ in range(3)]
