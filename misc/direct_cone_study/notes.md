# Combined direct/diffuse (forward-peak cone) investigation notes

## 2026-11-06 box cloud (ex_pprts_box_cld, 17x33x20, dx=200, dz=100, cloud tau~11 w0~0.9 g=0.9), SZA 70
- mcdmda (fixed) reference vs tenstream 3_10 (online bmc, tight tol): shadow total down mc 53 vs ts 85 (cone0) .. 99 (cone20).
- tenstream shadow edge is smeared over ~8 px (numerical diffusion of the 3-stream direct beam at SZA 70);
  mcdmda shadow edge is sharp (one px). This dominates the error (rmse total down ~64 W/m2 independent of cone).
- => forward-peak cone effect cannot be judged at SZA 70 with this setup; go to SZA 0 where the direct beam is not diffused.

## SZA 0, same box cloud (single cloud layer tau 10, 1 km clear air below)
run      shadow_tot  lit_tot  rmse_tot   (mcdmda: 168 / 505)
ds        225        489      30.8   (standard delta scaling)
cone0     202        493      21.2
cone5     203        492      21.6
cone10    206        492      22.8
cone20    216        489      27.1
- every cone makes it worse; delta scaling is worst.
- profile below cloud: mc 220 vs ts 275 at the centre, ts too low in the far field (465 vs 476):
  tenstream transports diffuse light in the clear air below the cloud diffusively (direction is reset at every voxel face),
  mcdmda ballistically (long grazing paths) -> ts keeps the transmitted light too concentrated below the cloud.
  This is the same 'directional memory loss' as the direct-beam smearing at SZA 70, but for the diffuse field.
- the cone does not address this; the box cloud has only one cloud layer, so no vertical in-cloud voxel boundaries
  where the forward peak could be kept.

## 1D slab check (ex_pprts_ex1, cloud one layer tau 10 w0 .995 g .9, SZA 0)
srfc edn: disort 889.6, twostr 911.8, mcdmda 903.7 (300k ph/px) -> mcdmda ok within 1.6%
next: multi-layer slab (10 layers tau 1 each) = clean test of the forward-peak cone (only vertical voxel boundaries)

## 1D reference check (tau 10, w0 .995, g .9, SZA 0, albedo .1)          R_toa  T_srfc  A_atm
independent numpy 1D MC (bmc_tools/mc1d.py, 4e5 photons)                30.82  66.51   9.32
mcdmda (single column of 3x3 homogeneous domain)                         31.06  66.30   9.27
DISORT (tenstream wrapper, 16..128 streams, converged)                   29.6   65.2   11.7   <-- absorption 25% too high!
twostream                                                                28.38  66.85  11.45
tenstream 3_10 online bmc (single cloud layer, dx 200 dz 100)            32.72  63.77   9.89
g=0 and w0=1 variants: DISORT and MC agree -> DISORT problem only with strong forward scattering + absorption
(possible issue in tenstream's disort wrapper / DISORT_rrtmg; follow up)
- CAUSE: disort/tenstream_disort.F90 set pmom(0)=1, pmom(1)=g, rest 0 (phase fct 1+3g cos, negative for g>1/3).
  FIXED (commit a6a2d881 on FJ/direct_cone): pmom(l)=g**l. Now DISORT 31.03/66.24/9.35 vs numpy MC 30.82/66.51/9.32; SZA60 49.93/44.37/10.13 vs 49.75/44.57/10.14.
- ml*_disort runs from run_ml.sh at SZA 0 were done with the OLD disort -> rerun.

## multi-layer slab (10 layers tau 1, g .9, w0 .995, dz 100, dx 200; tenstream 3_10 online bmc), % of S0*mu
SZA 0:  ref(disort fixed) R 30.99 T 66.30 | ds 35.71/61.39 | cone0 48.18/46.03 | 5: 47.05/47.34 | 10: 44.66/50.19 | 20: 39.69/56.25
SZA 60: ref 49.89/44.42 | ds 50.53/44.08 | cone0 55.46/38.75 | 20: 53.55/40.37
-> without forward-peak treatment the angular reset at each voxel face makes g=.9 behave nearly isotropic: R +17 pts!

## toy model (bmc_tools/toy.py): MC with tenstream-like angular reset at voxel faces
reproduces tenstream within ~1 pt (cone0 48.8/45.4, delta 36.1/61.1, cone20 40.4/55.3). 6 s per scan.
predictions (dR vs exact, SZA0 / SZA60): delta +5.2/+1.4; cone30 ph0 +5.3/+3.6; cone45 ph0 +1.2/-0.1; cone60 ph0 -0.9/-3.4;
cone30 ph1 (diffuse->direct only after scattering in box) +0.9/+2.1; ph1 incl. unscattered diffuse overcorrects strongly.

## tenstream confirms toy for wider cones (multi-layer slab, % of S0*mu)
SZA 0:  ref 30.99/66.30 | ds 35.71/61.39 maxerr 5.0 | cone30 35.57/61.36 maxerr 5.0 | cone45 31.30/66.27 maxerr 1.0 (!) | cone60 28.95/69.81 maxerr 3.5
SZA 60: ref 49.89/44.42 | ds 50.53/44.08 maxerr 3.9 | cone30 52.68/41.42 maxerr 5.0 | cone45 49.10/44.83 maxerr 8.5 (totals ok, but too deep penetration in upper cloud) | cone60 46.33/47.91 maxerr 11
- layer absorption error (toy, rms over cloud layers, % of S0*mu; layer abs ~1%):
  SZA0: delta .055 c30p0 .033 c45p0 .049 c30p1 .088 ; SZA60: delta .090 c30p0 .070 c45p0 .101 c30p1 .074
  -> heating-rate wise, all are close; c30p0 slightly best. R/T wise c45p0 best at SZA0.

## toy robustness (72 cases: g .7/.85/.9, tau_lay .2/1/4, SZA 0/30/60/75, dx 200/1000), mean max(|dR|,|dT|) pts
overall: cone0 5.24, delta 2.15, c30p0 2.12, c45p0 2.41, c60p0 3.48, c30p1 2.40; worst: 21.7/7.7/6.8/10.2/12.4/8.6
by SZA: cones win at SZA 0-30 (c45p0 1.63 vs delta 2.74 at SZA0), delta wins at SZA 60-75.
## tenstream: delta + cone (slab)
SZA0: ds 35.71/61.39 err 5.0 | ds+c10 35.39/61.64 err 4.8 | ds+c20 34.69/63.03 err 3.5
SZA60: ds 50.53/44.08 err 3.9 | ds+c10 50.48/43.94 err 3.8 | ds+c20 50.35/44.07 err 4.0   -> safe, modest gain

## large toy study (slurm, 120 cases, 40k photons; g .85/.9, tau_lay .1/.3/1/3, SZA 0/20/40/60/75, dx 100/200/1000) -> toy/summary.txt
RT = max(|dR|,|dT|) pts, AB = rel rms layer absorption error %; noise floor RT .2, AB 1.2
delta 2.70/8.16 | cone0 7.20/13.93 | c20p0 3.68/6.90 | c25p0 3.10/6.63 | c30p0 2.73/7.09 | c35p0 2.58/7.87
d+c20 2.36/8.80 | d+c30 2.09/9.86 | c20p1 2.80/6.50
-> realistic gain vs delta scaling ~20% (RT with d+c30, AB with c20p1 / c25p0). Cones help at SZA<=40, not at 60-75.
-> remaining AB errors 6-8% are due to the angular reset of the *diffuse* field, not the forward peak.

## other stream layouts (toy + tenstream)
- 3_24 == 3_10 in the slab (tenstream: identical within noise): its extra top streams are azimuthal quadrants,
  no polar resolution -> no help for the forward peak.
- 3_30 has a polar cap mu>=0.75 around each face normal + 4 quadrants. toy: in g=.9 w0=.995 slab absorption 16.7% vs 9.3% exact
  (photons in the side-face caps are re-emitted cos-weighted around the horizontal normal -> drift to horizontal paths, trapping).
  boxmc 3_30 emits the same way (mu = sqrt(U(alim^2,1))) -> expect real 3_30 to overestimate absorption too (checking).
- 'sun cap' diffuse stream (cap around the sun direction, resampled within the cap): worse than the collimated cone
  (SZA0: dR +9.1 vs cone25 +6.9; with delta +4.3/AB 7.7 vs delta 5.2/5.9) -> collimating into the direct stream is the better idea.

## 3D box cloud, vertically extended (10 layers tau 1, 2.2x2.2 km, 1 km thick, 500 m above srfc), mcdmda 4x30k ph/px
SZA 0  (mc shadow 163.3 lit 513.9): ds 217/489 rmse 34.0 | cone0 145/494 rmse 19.9 | cone20 186/493 21.8 | cone30 213/491 31.4 | ds+c20 229/488 39.0
SZA 30 (mc shadow 168.6 lit 521.0): ds 239/489 55.2 | cone0 168/497 51.3 | cone20 213/492 49.2 | cone30 240/489 55.2 | ds+c20 249/487 58.3
-> in 3D, collimation (delta scaling or wide cones) puts too much light straight through the cloud into the shadow
   (re-collimation error: real forward-scattered photons spread laterally by ~25 deg per scattering).
   1D wants the forward peak kept, 3D wants it not collimated -> trade-off; small cones (5-10 deg) being tested.
- tenstream 3_30 (azimuth 225 to avoid unimplemented dir3_to_diff30_coeff_symmetry): slab SZA0 ds: R 36.4 T 56.2 A 13.1 (ref 31.0/66.3/9.3)
  -> 3_30 overestimates absorption by 40% in g=.9 clouds, as predicted by the toy (side-face caps trap photons).
- tenstream 3_30 tight tol (slab, A ref 9.3/10.1): SZA0 ds A 14.6 cone0 16.5 cone30 14.6; SZA60 ds 16.2 cone0 15.2.
  toy predicted 3_30 cone0 R 45.1 T 42.4 vs tenstream 45.4/42.4 -> design property of the 3_30 stream layout
  (side-face cap trapping), not a code bug. 3_30 is unsuitable for strongly forward-scattering clouds.
- mcdmda robustness: aborts with 'waited 5 s but the queues havent freed up' on loaded nodes; use -mcdmda_queue_waittime 1800.

## 3D v2 absorption profiles, SZA 0 (ref: 3x mcdmda 8 ranks x 30k ph/px, run-to-run noise 0.04%)
run      dTOAdom dTOAcld | dAbs_dom dAbs_cld | rmsP_dom% rmsP_cld%
ds        +6.1   +39.9   |  +0.60   +3.27    |  8.4  10.1
cone0    +16.0   +71.9   |  +1.05   +5.80    | 23.6  26.3
cone5    +15.0   +68.5   |  +1.24   +5.88    | 23.5  24.6
cone10   +13.2   +61.9   |  +1.31   +5.70    | 21.5  21.2
cone20    +9.2   +48.6   |  +1.06   +5.13    | 15.7  16.8
cone30    +6.1   +38.2   |  +0.65   +4.12    |  9.2  12.8
dc20      +5.5   +37.5   |  +0.18   +1.43    |  3.2   4.9   <- best heating rates
=> conflict: surface shadow prefers narrow cone (5-10), cloud heating/TOA prefer delta+cone20.
   All variants overestimate cloud-top reflectance strongly (+38..+72 W/m2 on ref 81).
SZA30/60 mcdmda failed: UCX sysv shm exhausted (lateral photon traffic). Resubmitted with UCX_TLS=self,tcp (job 17806923).

## toy: truncation strength f=g^p (g .9, dx 200, 12 cases, 20k ph) RT[pts] / AB[%]
delta(p2) 2.98/8.24; p1.0 3.60/15.7; p1.25 3.16/12.6; p1.5 2.91/10.7; p1.5+c10 2.80/10.6;
p2+c20 2.70/8.9; p3+c20 2.89/6.9.  => stronger truncation trades R bias for worse layer absorption; no free lunch.

## 3D v2 absorption profiles SZA30/60 (mcdmda over tcp: 450 s instead of 5800 s!)
cloud-column rel rms heating-profile error [%]   SZA0   SZA30  SZA60
  delta                                           10.1    6.2    7.9
  cone0                                           26.3   20.1    7.4
  cone10                                          21.2   18.1    9.5
  cone20                                          16.8   13.3   10.4
  cone30                                          12.8   10.9    9.9
  delta+cone20                                     4.9    6.4    8.4
TOA eup over cloud overestimated by all (+38..+80 W/m2), worst at SZA60.
=> heating rates: delta+cone20 halves the error at overhead sun, neutral at SZA30/60. Cone without delta never wins.
   surface shadow (earlier): cone 5-10 at SZA0/30, cone 30 at SZA60. No single setting optimal for both.

## toy: face-angle-preserving diffuse layouts (2026-11-07), reset="rings" in toy.py
top/bot faces: nring polar rings around the normal (equal flux width in mu^2), photon keeps its ring.
side faces: elevation bands in |d_z| (equal width), keep sign of d_z and of normal comp (nband=1 == 3_10).
Rings around the side-face normal do NOT work (scramble mu_z): r32 with side rings R 38.6 vs exact 30.9.
Diffuse stream count = 2*nring + 8*nband.  12 slab cases g .9 dx 200, 20k ph:  RT[pts] / AB[%]
                    plain        +delta       +cone20
  r1e1 (3_10,n10)  7.69/15.4    2.96/8.50    4.12/7.28
  r2e2 (n20)       4.69/10.1    2.39/7.77    2.32/5.03
  r3e2 (n22)       3.58/8.15    2.08/7.55    1.63/4.11
  r4e2 (n24)       2.90/5.72    1.77/7.18    1.29/3.52
  r4e4 (n40)       1.81/3.64    1.75/7.35    0.69/2.91
  r8e8 (n80)       0.96/2.03    1.68/7.19    0.51/2.68
=> delta scaling saturates at AB ~7% (its own error); rings+cone keep converging.
   r4e2+cone20 (24 streams, same count as 3_24) halves 3_10+delta errors (1.29/3.52 vs 2.96/8.50).
   equal-angle spacing worse than equal-flux; r4e2 cone scan: c10 1.99/4.04, c15 1.44/3.38, c20 1.29/3.52, c25 1.22/3.95, c30 1.26/4.92.
   caveat: slab only (periodic); lateral 3D transport of elevation-band side streams untested.

## toy: what can be gained at 10 streams? (2026-11-08)
(a) stream shape: 3_10 sectors with non-lambertian intra-stream pdf (top/bot mu^p, sides x |d_z|^q), reset="shaped".
    Reduces R bias strongly (SZA0 R 48.8 -> 34.2 at p3q2, exact 30.9) but layer absorption gets worse in every combo:
    best shaped: p2q0+delta 2.74/11.5, p2q0+c20 2.90/9.2  vs  3_10+delta 2.98/8.2, 3_10+c20 4.22/7.6, 3_10+d+c20 2.70/8.9.
    A fixed shape is wrong somewhere in the profile (too much penetration).
(b) redistribute 10 streams (3 rings top/bot, unsplit sides): catastrophic (AB 129%), side up/down split is essential.
=> at 10 streams the layout/shape gives no benefit; the only levers are delta and the direct cone (already quantified).

## coupled direct/diffuse solve (phase 1) implemented, commit 2d009d84 on FJ/direct_cone (2026-11-09)
-pprts_diff2dir (+ -bmc_direct_cone_deg, -bmc_online), block Gauss-Seidel over solve_direct/solve_diffuse, 2-4 iterations.
Validation: slab SZA0 c20 coupled R/T/Tdir 36.87/59.68/18.71 vs toy phase1 37.10/59.57/18.23 (matches);
energy conservation w0=1 within +-0.1% (KSP precision); uncoupled path bit identical (mock LUT, old vs new binary);
azimuth convention control (flipped phi at SZA60) is worse -> convention right.
3D box cloud (box3d_cpl) vs mcdmda:
  cloud heating rms%  SZA0  SZA30 SZA60 | surface rmse tot   SZA0  SZA30 SZA60
  delta               10.1   6.2   7.9  |                    33.7  55.2  63.8
  cone20               16.8  13.3  10.4 |                    22.0  49.2  64.9
  cpl20                13.4   9.9   8.4 |                    30.3  53.3  59.3
  delta+cone20 (dc20)   4.9   6.4   8.4 |                    38.9  58.3  65.0
  delta+cpl20           8.5   6.0   7.6 |                    41.4  59.9  64.9
  cone10               21.2  18.1   9.5 |                    18.1  48.0  79.5
  cpl10                20.0  15.9   8.5 |                    18.9  47.9  76.8
=> coupling improves the cone for heating at all SZA (20-35%), improves surface only at SZA60;
   at SZA0/30 coupling brightens the shadow (converted diffuse goes straight down instead of spreading).
   No coupled variant clearly beats delta (+cone20 phase0) overall; best-of-class only marginal (dcpl20 at SZA30/60 heating).
