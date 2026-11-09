# Direct-cone / diffuse-to-direct coupling study (archived, not for merge)

Research on keeping the forward scattering peak by counting radiation within a cone
around the sun as direct radiation (`-bmc_direct_cone_deg`, phase 0) and coupling
diffuse radiation scattered into the cone back into the direct streams
(`-pprts_diff2dir`, phase 1, block Gauss-Seidel between direct and diffuse solve).

- `notes.md`: chronological study notes and results (vs. mcdmda reference)
- `tools/`: python toy model (`toy.py`, MC with tenstream-like face reset) and scan / analysis scripts
- `slurm/`: job scripts and run lists for the slab, 3D box cloud and toy experiments

Conclusion: in 3D neither the cone nor the coupled solve clearly beats delta scaling
(+ cone 20 deg phase 0) for the 3_10 solver. Branch kept for archival only.
