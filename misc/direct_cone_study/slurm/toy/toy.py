"""MonteCarlo emulation of the voxel-face 'angular reset' of a tenstream-like stream discretization
in a horizontally homogeneous (periodic) layered atmosphere.

photons move through voxels (dx, dy, dz[k]); at every voxel face they cross:
  - direct photons keep the solar direction,
  - with a cone: photons whose direction is within the cone half angle of the sun (and that exit through a face
    that the direct beam can reach, i.e. not upwind) become direct (dir := sun)
      phase0: only photons of direct lineage (were direct when they entered the voxel)
      phase1: also diffuse photons ('diff2dir'), optionally only if they scattered inside this voxel
  - diffuse photons get a new direction, cosine-weighted around the face normal, restricted to the
    sector of their stream (3_10-like: top/bot faces one sector per hemisphere, side faces split in
    upward/downward going halves)
with reset='none' this is an exact 3D MonteCarlo of the plane-parallel problem (reference).

usage as module: run(cfg) -> dict(R, T, Tdir, A)
"""
import numpy as np


def hg(r, g):
    if abs(g) < 1e-6:
        return 2 * r - 1
    return (1 + g * g - ((1 - g * g) / (1 - g + 2 * g * r)) ** 2) / (2 * g)


def scatter(d, g, rng):
    n = d.shape[0]
    ct = hg(rng.random(n), g)
    st = np.sqrt(np.maximum(0, 1 - ct * ct))
    phi = 2 * np.pi * rng.random(n)
    # orthonormal basis around d
    a = np.where(np.abs(d[:, 2:3]) < 0.9, np.array([[0, 0, 1.0]]), np.array([[1.0, 0, 0]]))
    u = np.cross(d, a)
    u /= np.linalg.norm(u, axis=1)[:, None]
    v = np.cross(d, u)
    return ct[:, None] * d + st[:, None] * (np.cos(phi)[:, None] * u + np.sin(phi)[:, None] * v)


def lambert(n, rng):
    """cosine weighted directions around +z"""
    mu = np.sqrt(rng.random(n))
    phi = 2 * np.pi * rng.random(n)
    s = np.sqrt(1 - mu * mu)
    return np.stack([s * np.cos(phi), s * np.sin(phi), mu], axis=1)


def reset_dir(d, axis, rng, split_sides=True):
    """resample direction cos-weighted around the face normal (axis), keep sign along the axis,
    for side faces keep the sign of the vertical component (3_10-like up/down split of side streams)"""
    n = d.shape[0]
    out = np.empty_like(d)
    for _ in range(100):
        l = lambert(n, rng)
        # map local frame: normal = axis
        if axis == 2:
            new = l.copy()
        elif axis == 0:
            new = np.stack([l[:, 2], l[:, 0], l[:, 1]], axis=1)
        else:
            new = np.stack([l[:, 0], l[:, 2], l[:, 1]], axis=1)
        new[:, axis] *= np.sign(d[:, axis])
        if axis != 2 and split_sides:
            new[:, 2] = np.abs(new[:, 2]) * np.where(d[:, 2] >= 0, 1, -1)
        out = new
        break
    return out


def run(tau_lay, w0_lay, g_lay, dz, dx, sza, albedo, N=100000, reset="sector", cone_deg=0.0,
        phase1=False, phase1_scattered_only=True, seed=1, phi0=0.0):
    rng = np.random.default_rng(seed)
    tau_lay, w0_lay, g_lay, dz = map(np.asarray, (tau_lay, w0_lay, g_lay, dz))
    nlay = len(tau_lay)
    kext = tau_lay / dz
    th, ph = np.deg2rad(sza), np.deg2rad(phi0)
    sun = np.array([np.sin(th) * np.sin(ph), np.sin(th) * np.cos(ph), -np.cos(th)])  # direction of travel, z up
    cos_cone = np.cos(np.deg2rad(cone_deg))
    L = np.array([dx, dx, 0.0])

    pos = np.zeros((N, 3))
    pos[:, 0] = rng.random(N) * dx
    pos[:, 1] = rng.random(N) * dx
    lay = np.zeros(N, int)  # layer index from top
    pos[:, 2] = dz[0]  # local z within layer, 0 = bottom, dz = top
    d = np.tile(sun, (N, 1))
    w = np.ones(N)
    direct = np.ones(N, bool)
    lin_dir = np.ones(N, bool)  # direct lineage inside the current voxel (was direct at voxel entry)
    scat_here = np.zeros(N, bool)
    alive = np.ones(N, bool)
    R = T = Tdir = A = 0.0
    Alay = np.zeros(nlay)

    tau_left = -np.log(rng.random(N))
    it = 0
    while alive.any():
        it += 1
        idx = np.nonzero(alive)[0]
        p, dd, k = pos[idx], d[idx], lay[idx]
        # distance to faces
        with np.errstate(divide="ignore", invalid="ignore"):
            tx = np.where(dd[:, 0] > 0, (dx - p[:, 0]) / dd[:, 0], np.where(dd[:, 0] < 0, -p[:, 0] / dd[:, 0], np.inf))
            ty = np.where(dd[:, 1] > 0, (dx - p[:, 1]) / dd[:, 1], np.where(dd[:, 1] < 0, -p[:, 1] / dd[:, 1], np.inf))
            tz = np.where(dd[:, 2] > 0, (dz[k] - p[:, 2]) / dd[:, 2], np.where(dd[:, 2] < 0, -p[:, 2] / dd[:, 2], np.inf))
        texit = np.minimum(np.minimum(tx, ty), tz)
        s = tau_left[idx] / kext[k]
        sc = s < texit
        # --- scattering events
        i_s = idx[sc]
        if i_s.size:
            pos[i_s] += d[i_s] * s[sc][:, None]
            kk = lay[i_s]
            A += (w[i_s] * (1 - w0_lay[kk])).sum()
            np.add.at(Alay, kk, w[i_s] * (1 - w0_lay[kk]))
            w[i_s] *= w0_lay[kk]
            # group by g value
            for gv in np.unique(g_lay[kk]):
                m = g_lay[kk] == gv
                d[i_s[m]] = scatter(d[i_s[m]], gv, rng)
            direct[i_s] = False
            scat_here[i_s] = True
            tau_left[i_s] = -np.log(rng.random(i_s.size))
        # --- face crossings
        i_f = idx[~sc]
        if i_f.size:
            t = texit[~sc]
            tau_left[i_f] -= t * kext[lay[i_f]]
            pos[i_f] += d[i_f] * t[:, None]
            ax = np.argmin(np.stack([tx[~sc], ty[~sc], tz[~sc]], axis=1), axis=1)
            # cone refill / reset of diffuse photons
            if reset != "none":
                cosang = d[i_f] @ sun
                downwind = d[i_f, ax] * sun[ax] > 0
                incone = (cosang > cos_cone) & downwind & (~direct[i_f]) & (cone_deg > 0)
                if phase1:
                    ok = lin_dir[i_f] | (scat_here[i_f] if phase1_scattered_only else True)
                else:
                    ok = lin_dir[i_f]
                ref = incone & ok
                d[i_f[ref]] = sun
                direct[i_f[ref]] = True
                rs = ~direct[i_f]
                for a in range(3):
                    m = rs & (ax == a)
                    if m.any():
                        d[i_f[m]] = reset_dir(d[i_f[m]], a, rng)
            # move to neighbour
            for a in (0, 1):
                m = ax == a
                ii = i_f[m]
                pos[ii, a] = np.where(d[ii, a] > 0, 0.0, dx)
            m = ax == 2
            ii = i_f[m]
            up = d[ii, 2] > 0
            # leaving through top
            iu = ii[up]
            top_toa = lay[iu] == 0
            R += w[iu[top_toa]].sum()
            alive[iu[top_toa]] = False
            iu = iu[~top_toa]
            lay[iu] -= 1
            pos[iu, 2] = 0.0
            # leaving through bottom
            idn = ii[~up]
            srf = lay[idn] == nlay - 1
            isf = idn[srf]
            T += w[isf].sum()
            Tdir += w[isf[direct[isf]]].sum()
            w[isf] *= albedo
            d[isf] = lambert(isf.size, rng)
            direct[isf] = False
            pos[isf, 2] = 0.0  # stays in lowest layer, at its bottom, going up
            idn = idn[~srf]
            lay[idn] += 1
            pos[idn, 2] = dz[lay[idn]]
            # entering a new voxel: update lineage
            lin_dir[i_f] = direct[i_f]
            scat_here[i_f] = False
        kill = alive & (w < 1e-6)
        alive &= ~kill
        if it > 200000:
            raise RuntimeError("too many iterations")
    return dict(R=R / N, T=T / N, Tdir=Tdir / N, A=A / N, Alay=Alay / N)


def slab(ncld=10, tau=10.0, w0=0.995, g=0.9, nclr_top=5, nclr_bot=5, dz=100.0, delta=False, **kw):
    taus = [0.005 / 20] * nclr_top + [tau / ncld] * ncld + [0.005 / 20] * nclr_bot
    w0s = [1.0] * nclr_top + [w0] * ncld + [1.0] * nclr_bot
    gs = [0.0] * nclr_top + [g] * ncld + [0.0] * nclr_bot
    taus, w0s, gs = map(np.array, (taus, w0s, gs))
    if delta:
        f = gs**2
        taus = taus * (1 - w0s * f)
        w0s = w0s * (1 - f) / (1 - w0s * f)
        gs = (gs - f) / (1 - f)
    return run(taus, w0s, gs, np.full(len(taus), dz), **kw)


if __name__ == "__main__":
    import sys

    sza = float(sys.argv[1]) if len(sys.argv) > 1 else 0.0
    N = int(float(sys.argv[2])) if len(sys.argv) > 2 else 50000
    base = dict(dx=200.0, sza=sza, albedo=0.1, N=N)
    print("exact (no reset)", slab(reset="none", **base))
    print("reset cone0     ", slab(**base))
    print("reset delta     ", slab(delta=True, **base))
    for c in (5, 10, 20, 30, 45):
        print(f"reset cone{c:<3d} ph0", slab(cone_deg=c, **base))
