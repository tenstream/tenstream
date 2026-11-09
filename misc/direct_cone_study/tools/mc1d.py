"""independent 1D plane-parallel MonteCarlo (vectorized numpy), HG phase function, lambertian surface
single homogeneous layer of optical thickness tau, ssa w0, asymmetry g, sun at zenith angle sza.
returns reflectance at TOA, transmittance (total down at the surface) and atmospheric absorptance.
photon weights: absorption by survival weight (w0 per scattering event), surface reflection by albedo weight.
"""
import sys
import numpy as np


def hg(r, g):
    if abs(g) < 1e-6:
        return 2 * r - 1
    return (1 + g * g - ((1 - g * g) / (1 - g + 2 * g * r)) ** 2) / (2 * g)


def run(tau, w0, g, sza, albedo, N, seed=1):
    rng = np.random.default_rng(seed)
    mu = np.full(N, np.cos(np.deg2rad(sza)))  # direction cosine, positive = downward
    z = np.zeros(N)  # optical depth from top
    w = np.ones(N)
    R = T = A = Asrfc = 0.0
    alive = np.ones(N, bool)
    while alive.any():
        idx = np.nonzero(alive)[0]
        s = -np.log(rng.random(idx.size))
        znew = z[idx] + s * mu[idx]
        top = znew < 0
        bot = znew > tau
        R += w[idx[top]].sum()
        alive[idx[top]] = False
        # surface hit: tally transmitted, reflect with albedo
        ib = idx[bot]
        T += w[ib].sum()
        Asrfc += (w[ib] * (1 - albedo)).sum()
        w[ib] *= albedo
        z[ib] = tau
        mu[ib] = -np.sqrt(rng.random(ib.size))  # lambertian upward
        # scattering inside
        isc = idx[~top & ~bot]
        z[isc] = znew[~top & ~bot]
        A += (w[isc] * (1 - w0)).sum()
        w[isc] *= w0
        cost = hg(rng.random(isc.size), g)
        phi = 2 * np.pi * rng.random(isc.size)
        m = mu[isc]
        mu[isc] = m * cost + np.sqrt(np.maximum(0, 1 - m * m)) * np.sqrt(np.maximum(0, 1 - cost * cost)) * np.cos(phi)
        # kill negligible weights (roulette free: just stop below 1e-9)
        dead = w < 1e-9
        A += w[dead & alive].sum() * 0  # negligible
        alive &= ~dead
    return R / N, T / N, A / N, Asrfc / N


if __name__ == "__main__":
    tau, w0, g, sza, alb = map(float, sys.argv[1:6])
    N = int(float(sys.argv[6])) if len(sys.argv) > 6 else 200000
    R, T, A, As = run(tau, w0, g, sza, alb, N)
    print(f"tau {tau} w0 {w0} g {g} sza {sza} alb {alb} N {N}: R {100*R:.2f}  T(total down srfc) {100*T:.2f}  A_atm {100*A:.2f}  A_srfc {100*As:.2f}  sum {100*(R+A+As):.2f}")
