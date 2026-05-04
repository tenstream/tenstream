# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What TenStream Is

TenStream is a parallel 3D radiative transfer solver for atmospheric heating rates. It solves the radiative transfer equation over structured (DMDA) and unstructured (PETSc PLEX) grids using MPI, and is coupled to NWP and LES models (COSMO, DALES, PALM, UCLA-LES, WRF). The solver name comes from the default configuration of 3 direct + 10 diffuse streams.

## Build

Out-of-source builds are **required**. Compiler is detected from the `FC` environment variable to auto-select a config from `config/<syst>.cmake`.

```bash
export FC=mpif90 CC=mpicc CXX=mpicxx
export PETSC_DIR=/path/to/petsc PETSC_ARCH=default

mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=RELEASE ..
make -j$(nproc)
```

To build dependencies (PETSc + NetCDF) from scratch:
```bash
FC=mpif90 CC=mpicc CXX=mpicxx misc/build_dependencies.sh
```

Key CMake options:
- `-DPETSC_INSTALL=<path>` — path to PETSc if not using env vars
- `-DENABLE_RRTM=ON/OFF`, `-DENABLE_PLEXRT=ON/OFF`, `-DENABLE_ECCKD=ON/OFF`, `-DENABLE_SPECINT=ON/OFF`, `-DENABLE_REPWVL=ON/OFF` — toggle optional modules
- `-DENABLE_ALL=ON` — enable all optional modules
- `-DCMAKE_BUILD_TYPE=DEBUG` — enables bounds checking, backtraces
- `-DENABLE_GCOV=ON` — code coverage flags
- `SYST=gcc|intel|lrz|mistral|...` — load a system-specific config from `config/`

## Testing

Tests use the **pFUnit** framework (MPI-enabled), managed as a CMake external project under `external/pfunit`.

```bash
# Run all tests except slow boxmc/LUT tests (default CI config):
cd build && make check

# Run a specific test:
ctest --output-on-failure -R test_pprts_symmetry

# Run LUT tests only:
make checkLUT

# Run all tests verbosely:
make check_verbose
```

Each test lives in `tests/<test_name>/` with pFUnit annotations (`@test`, `@before`, `@after`, `@mpiTest`). The `NUMPROC` variable in `test_settings.cmake` controls MPI rank count.

LUT (look-up table) tests require downloading binary files first:
```bash
make download_common_luts   # downloads ~3 most-used LUTs
make download_all_luts      # downloads everything
```

## Code Formatting

All Fortran source must be formatted with **fprettify** before committing. CI enforces this.

```bash
# Format all source files in-place:
maint/fprettify.sh

# Check (diff only, no changes):
maint/fprettify.sh -d

# Format only staged files:
maint/fprettify.sh -s
```

Settings: 2-space indent, 192-char line length, whitespace level 3, strict indent, lowercase keywords.

## Architecture

### Module naming convention
All Fortran modules are prefixed `m_`. Example: `m_pprts`, `m_optprop`, `m_helper_functions`.

### Source layout
| Directory | Purpose |
|-----------|---------|
| `src/` | Core solver — boxmc ray tracer, pprts solver, optical properties, helpers |
| `plexrt/` | RT on unstructured grids (ICON/DALES) via PETSc PLEX |
| `rrtmg/` | RRTMG gas optics (LW + SW, includes upstream Fortran source) |
| `ecckd/` | ECCKD gas optics |
| `specint/` | Spectral integration framework |
| `repwvl/` | Representative wavelength approach |
| `disort/` | Coupled DISORT 1D solver |
| `c_wrapper/` | C API (`f2c_pprts`) for calling Fortran solver from C |
| `tests/` | pFUnit test suites, one subdirectory per suite |
| `examples/` | Runnable examples; `ex_*.F90` contain `main()`, `*.F90` are library code |
| `config/` | Compiler/system-specific CMake flags (gcc, intel, lrz, mistral, …) |
| `external/` | fypp preprocessor (submodule), pFUnit (CMake external project) |
| `misc/` | Build helper scripts, LUT download script, benchmarks |

### Solver stream configurations
Solvers are typed as `t_solver_N_M` where N = direct streams, M = diffuse streams. Available types (defined in `m_pprts_base`): `t_solver_1_2`, `t_solver_3_6`, `t_solver_3_10`, `t_solver_3_16`, `t_solver_3_24`, `t_solver_3_30`, `t_solver_8_10`, `t_solver_8_16`, `t_solver_8_18`. Wedge variants exist for `plexrt`. Higher stream counts give better accuracy at higher cost.

### Core data flow
1. **Optical properties** (`m_optprop`, `m_optprop_LUT`, `m_optprop_ANN`): transfer coefficients for each grid cell come from Look-Up Tables (`.mmap4` binary files) built offline by `src/createLUT.F90` using the BoxMC Monte Carlo ray tracer (`m_boxmc`). An ANN alternative exists for `3_10`.
2. **Solver** (`m_pprts`, `m_pprts_base`): assembles a sparse linear system via PETSc DMDA, sets up KSP solvers (separate for direct `dir_` and diffuse `diff_` components), and solves for irradiance.
3. **Gas optics** (`rrtmg/`, `ecckd/`, `specint/`): compute per-band optical properties from atmospheric state, then call the core solver for each band/g-point.

### fypp preprocessing
Files ending in `.fypp` are Fortran templates processed by [fypp](https://github.com/aradi/fypp) at build time into `.F90`. Used for rank-generic array utilities (`helper_functions.fypp`, `netcdfio.fypp`, `search.fypp`, `intersection.fypp`, `linked_list.fypp`). The preprocessor lives at `external/fypp`.

### Runtime configuration
All solver options are passed via the **PETSc options database** (command line or `~/.petscrc`). Key TenStream options (see `src/tenstream_options.F90`):

```
-lut_basename <path_prefix>      path to LUT binary files
-schwarzschild                   use Schwarzschild for thermal instead of twostream
-twostr_ratio <N>                use 1D twostream when aspect ratio dz/dx > N (default 2)
-calc_nca                        enable Non-local Convergence Adjustment
-skip_thermal                    zero out thermal flux/absorption
-topography                      enable surface topography via ray bending
-max_solution_err <W>            skip solve if estimated error below threshold
-show_options                    print all tenstream options
```

PETSc KSP solver options use prefixes `-dir_ksp_*` (direct) and `-diff_ksp_*` (diffuse). See `doc/petscrc_example` for tuning examples including AMG preconditioners.

### Key types in `m_pprts_base`
- `t_atmosphere` — atmospheric state (optical thickness, single scatter albedo, asymmetry factor per cell)
- `t_solver` — base solver type; extended by stream-specific types
- `t_state_container` — flux solution container
- `t_suninfo` — solar geometry (zenith/azimuth)
- `t_coord` — local/global DMDA index ranges and ghost-point layout

### Data types (`m_data_parameters`)
- `ireals` — working precision real (single or double depending on PETSc build)
- `irealLUT` — LUT storage precision
- `iintegers` — working precision integer (32 or 64 bit depending on PETSc build)
- `mpiint` — MPI integer type

## Supported Compilers
GFortran, Intel (`ifort`/`ifx`), XLF IBM. The `FC` environment variable selects the config. Intel oneAPI uses `mpiifx`/`mpiicx`.
