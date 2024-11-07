def eddington_ec(dtau, w0, g, mu0):
    import numpy as np
    dtau, w0, g, mu = [ np.atleast_1d(_) for _ in (dtau, w0, g, mu0) ]

    f = 0.75 * g
    g1 = 2. - w0 * (1.25 + f)
    g2 = w0 * (0.75 - f)
    g3 = 0.5 - mu0 * f

    dtau_slant = np.maximum(dtau / np.maximum(1e-30, mu0), 0)

    g4 = 1. - g3
    alpha1 = g1 * g4 + g2 * g3
    alpha2 = g1 * g3 + g2 * g4

    A = np.sqrt(np.maximum((g1 - g2) * (g1 + g2), 1e-12))
    k_mu0 = A * mu0
    k_g3 = A * g3
    k_g4 = A * g4

    e0 = np.exp(-dtau_slant)
    tdir = e0

    e = np.exp(-A * dtau)
    e2 = e * e
    k_2_e = 2 * A * e
    k_mu0[(k_mu0 > 1.-1e-8) & (k_mu0 < 1.+1e-8)] = 1. - 1e-7
    beta = 1. / (A + g1 + (A - g1) * e2)
    r = g2 * (1. - e2) * beta
    t = k_2_e * beta

    beta = w0 * beta / (1. - k_mu0 * k_mu0)

    sdir = beta * (k_2_e * (g4 + alpha1 * mu0) - e0 * ((1 + k_mu0) * (alpha1 + k_g4) - (1 - k_mu0) * (alpha1 - k_g4) * e2))
    rdir = beta * ((1 - k_mu0) * (alpha2 + k_g3) - (1 + k_mu0) * (alpha2 - k_g3) * e2 - k_2_e * (g3 - alpha2 * mu0) * e0)

    single_scat_mask = dtau_slant < 1e-6
    t   [single_scat_mask] = (1. - g1 * dtau         )[single_scat_mask]
    r   [single_scat_mask] = (g2 * dtau              )[single_scat_mask]
    sdir[single_scat_mask] = ((1. - g3) * w0 * dtau  )[single_scat_mask]
    rdir[single_scat_mask] = (g3 * w0 * dtau         )[single_scat_mask]
    tdir[single_scat_mask] = (1. - dtau_slant        )[single_scat_mask]

    return np.vstack([t, r, rdir, sdir, tdir]).T
