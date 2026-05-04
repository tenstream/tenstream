#!/usr/bin/env python3
"""Remove specific symbols from Fortran 'use module, only:' statements.

Usage: python3 remove_unused_imports.py
Applies the hardcoded removals map to source files in-place.
"""

import re
import sys
from pathlib import Path

# Map of filepath -> {module_lower -> [symbols_to_remove_lower]}
# All keys/values are lowercase for case-insensitive matching.
REMOVALS = {
    "disort/tenstream_disort.F90": {"m_data_parameters": ["mpiint"]},
    "ecckd/ecckd_optprop.F90": {
        "m_data_parameters": ["avogadro", "k_boltzmann", "r_dry_air"],
        "m_fu_ice": ["fu_ice_optprop", "fu_ice_data_solar", "fu_ice_data_thermal"],
    },
    "ecckd/ecckd_pprts.F90": {
        "m_data_parameters": ["default_str_len"],
        "m_pprts": ["pprts_get_result_tozero"],
    },
    "examples/libRadtran_cld_file/uvspec_cld_file.F90": {
        "m_dyn_atm_to_rrtmg": ["print_tenstr_atm"],
    },
    "examples/pprts/ex_pprts_specint_lw_sw.F90": {"m_netcdfio": ["ncwrite"]},
    "examples/pprts/ex_pprts_specint_tree.F90": {"m_helper_functions": ["tostr"]},
    "examples/pprts/pprts_box_cld.F90": {"m_tenstream_options": ["read_commandline_options"]},
    "examples/pprts/pprts_ex1.F90": {
        "m_data_parameters": ["zero", "pi"],
        "m_tenstream_options": ["read_commandline_options"],
    },
    "examples/pprts/pprts_repwvl_lw_sw.F90": {"m_data_parameters": ["zero", "one", "default_str_len"]},
    "examples/pprts/pprts_rrtm_iterations.F90": {"m_data_parameters": ["zero", "one"]},
    "examples/pprts/pprts_rrtm_lw_sw.F90": {"m_data_parameters": ["one", "default_str_len"]},
    "examples/pprts/pprts_specint_buildings.F90": {"m_data_parameters": ["default_str_len"]},
    "examples/pprts/pprts_specint_lw_sw.F90": {"m_data_parameters": ["zero", "one"]},
    "examples/pprts/pprts_specint_lw_sw_from_dump.F90": {
        "m_data_parameters": ["zero", "one", "pi"],
        "m_helper_functions": ["linspace", "spherical_2_cartesian"],
        "m_pprts": ["gather_all_tozero"],
        "m_dyn_atm_to_rrtmg": ["setup_tenstr_atm"],
        "m_petsc_helpers": ["getvecpointer", "restorevecpointer"],
    },
    "examples/pprts/pprts_specint_tree.F90": {
        "m_data_parameters": ["default_str_len"],
        "m_helper_functions": ["cstr"],
    },
    "examples/pprts_hill/pprts_hill.F90": {
        "m_data_parameters": ["zero", "one"],
        "m_helper_functions": ["itoa"],
    },
    "examples/wetterstein_ts/wetterstein.F90": {
        "m_tenstream_options": ["read_commandline_options"],
        "m_pprts": ["gather_all_tozero"],
    },
    "plexrt/icon_grid.F90": {"m_helper_functions": ["get_arg"]},
    "plexrt/icon_plex_utils.F90": {"m_data_parameters": ["one"]},
    "plexrt/nca_multi_tri.F90": {"m_helper_functions": ["chkwarn"]},
    "plexrt/plex_grid.F90": {
        "m_helper_functions": ["chkwarn", "rotation_matrix_local_basis_to_world", "strf2c"],
        "m_data_parameters": ["ireallut", "i7", "i8"],
        "m_icon_grid": ["iconull"],
    },
    "plexrt/plexrt_external_solvers.F90": {
        "m_helper_functions": ["chkwarn", "cstr"],
        "m_schwarzschild": ["b_eff"],
        "m_icon_plex_utils": ["dump_ownership", "dmplex_2d_to_3d"],
    },
    "plexrt/pprts2plex.F90": {"m_data_parameters": ["i0"]},
    "repwvl/fu_ice.F90": {"m_repwvl_base": ["t_repwvl_data"]},
    "repwvl/rayleigh.F90": {"m_data_parameters": ["iintegers", "ireal_dp"]},
    "repwvl/repwvl_plexrt.F90": {
        "m_data_parameters": ["default_str_len"],
        "m_helper_functions": ["get_arg"],
        "m_tenstream_options": ["read_commandline_options"],
        "m_mie_tables": ["destroy_mie_table"],
    },
    "repwvl/repwvl_pprts.F90": {
        "m_data_parameters": ["one", "default_str_len"],
        "m_pprts": ["pprts_get_result_tozero"],
    },
    "rrtmg/rrtmg/dyn_atm_to_rrtmg.F90": {"m_data_parameters": ["default_str_len"]},
    "rrtmg/rrtmg/pprts_rrtmg.F90": {
        "m_data_parameters": ["init_mpi_data_parameters", "i2", "i9", "pi"],
        "m_pprts": ["pprts_get_result_tozero"],
        "m_search": ["find_real_location"],
        "m_netcdfio": ["ncwrite"],
    },
    "src/buildings.F90": {
        "m_data_parameters": ["default_str_len"],
        "m_helper_functions": ["chkwarn", "get_arg"],
        "m_pprts_base": ["t_solver"],
    },
    "src/createLUT.F90": {"m_helper_functions": ["chkerr"]},
    "src/eddington.F90": {
        "m_data_parameters": ["iintegers", "pi"],
        "m_helper_functions": ["delta_scale_optprop"],
    },
    "src/interpolation.F90": {"m_data_parameters": ["mpiint"]},
    "src/mmap.F90": {"m_c_syscall_wrappers": ["c_lockf"]},
    "src/nca.F90": {"m_data_parameters": ["iintegers", "zero"]},
    "src/optprop_base.F90": {
        "m_helper_functions": ["ind_nd_to_1d"],
        "m_optprop_parameters": [
            "preset_param_phi11", "preset_aspect5", "preset_aspect11",
            "preset_aspect17", "preset_aspect90", "preset_w050", "preset_tau50",
        ],
    },
    "src/optprop.F90": {
        "m_helper_functions": ["deg2rad", "swap"],
        "m_data_parameters": ["i0", "i1", "inil"],
    },
    "src/optprop_LUT.F90": {
        "m_helper_functions": ["ind_1d_to_nd", "linspace", "mpi_logical_and", "ndarray_offsets", "rel_approx"],
        "m_data_parameters": ["i0", "i2", "i3", "i10", "nil", "inil", "imp_ireals"],
        "m_optprop_parameters": ["stddev_rtol"],
        "m_boxmc": ["t_boxmc"],
        "m_tenstream_interpolation": ["interp_4d"],
        "m_optprop_base": ["find_op_dim_by_name"],
    },
    "src/petsc_helpers.F90": {"m_data_parameters": ["i2", "i3"]},
    "src/pprts_base.F90": {
        "m_data_parameters": ["pi"],
        "m_helper_functions": ["imp_allreduce_min"],
    },
    "src/pprts_explicit.F90": {"m_helper_functions": ["mpi_logical_and"]},
    "src/pprts_external_solvers.F90": {
        "m_data_parameters": ["i4"],
        "m_helper_functions": ["chkwarn", "meanval"],
        "m_pprts_base": ["set_dmda_cell_coordinates"],
        "m_plex_grid": ["print_dmplex"],
    },
    "src/pprts.F90": {
        "m_data_parameters": ["imp_ireals", "i7", "i10"],
        "m_helper_functions": [
            "angle_between_two_normed_vec", "cross_3d", "deallocate_allocatable",
            "imp_allreduce_mean", "normalize_vec",
            "rotation_matrix_world_to_local_basis", "vec_proj_on_plane",
        ],
        "m_optprop": ["t_optprop", "t_optprop_cube"],
        "m_pprts_base": ["t_dof", "t_solver_log_events"],
    },
    "src/schwarzschild.F90": {"m_data_parameters": ["exp_minval"]},
    "src/tenstream_options.F90": {"m_data_parameters": ["ireallut"]},
    "src/twostream.F90": {"m_helper_functions": ["delta_scale_optprop"]},
    "src/xdmf_export.F90": {"m_helper_functions": ["ind_nd_to_1d", "ndarray_offsets"]},
    "specint/specint_plexrt.F90": {
        "m_helper_functions": ["get_arg", "tostr"],
        "m_data_parameters": ["iintegers"],
    },
    "tests/eddington/test_delta_eddington.F90": {"m_eddington": ["eddington_coeff_bm"]},
    "tests/test_boxmc_1_2/test_boxmc_1_2.F90": {
        "m_boxmc": ["t_boxmc"],
        "m_data_parameters": ["i1"],
        "m_helper_functions": ["tostr", "cstr"],
        "m_optprop_parameters": ["stddev_atol"],
        "m_eddington": ["eddington_coeff_zdun", "eddington_coeff_bm", "eddington_coeff_ec"],
    },
    "tests/test_boxmc_3_10/test_boxmc_3_10.F90": {
        "m_boxmc": ["t_boxmc"],
        "m_data_parameters": ["i1"],
        "m_optprop_parameters": ["stddev_atol"],
    },
    "tests/test_boxmc_3_24/test_boxmc_3_24.F90": {
        "m_boxmc": ["t_boxmc"],
        "m_data_parameters": ["i1"],
        "m_helper_functions": ["tostr", "deg2rad"],
        "m_optprop_parameters": ["stddev_atol"],
    },
    "tests/test_boxmc_3_30/test_boxmc_3_30.F90": {
        "m_boxmc": ["t_boxmc"],
        "m_data_parameters": ["i1"],
        "m_helper_functions": ["tostr", "deg2rad"],
        "m_optprop_parameters": ["stddev_atol"],
    },
    "tests/test_boxmc_3_6/test_boxmc_3_6.F90": {
        "m_boxmc": ["t_boxmc"],
        "m_data_parameters": ["i1"],
        "m_optprop_parameters": ["stddev_atol"],
    },
    "tests/test_boxmc_3_6_tau_scaling/test_boxmc_3_6_tau_scaling.F90": {
        "m_boxmc": ["t_boxmc"],
        "m_data_parameters": ["i1"],
        "m_optprop_parameters": ["stddev_atol"],
        "m_helper_functions": ["itoa"],
    },
    "tests/test_boxmc_8_10/test_boxmc_8_10.F90": {
        "m_boxmc": ["t_boxmc"],
        "m_data_parameters": ["i1"],
        "m_optprop_parameters": ["stddev_atol"],
    },
    "tests/test_boxmc_8_12/test_boxmc_8_12.F90": {
        "m_boxmc": ["t_boxmc"],
        "m_data_parameters": ["i1"],
        "m_optprop_parameters": ["stddev_atol"],
        "m_helper_functions": ["deg2rad"],
    },
    "tests/test_boxmc_8_16/test_boxmc_8_16.F90": {
        "m_boxmc": ["t_boxmc"],
        "m_data_parameters": ["i1"],
        "m_optprop_parameters": ["stddev_atol"],
        "m_helper_functions": ["deg2rad"],
    },
    "tests/test_boxmc_8_18/test_boxmc_8_18.F90": {
        "m_boxmc": ["t_boxmc"],
        "m_data_parameters": ["i1"],
        "m_optprop_parameters": ["stddev_atol"],
        "m_helper_functions": ["deg2rad"],
    },
    "tests/test_buildings/test_buildings.F90": {
        "iso_fortran_env": ["real32", "real64"],
        "m_data_parameters": ["default_str_len"],
    },
    "tests/test_convolution/test_convolution.F90": {"m_helper_functions": ["chkerr"]},
    "tests/test_geometric_coeffs/test_geometric_coeffs.F90": {
        "m_boxmc": ["t_boxmc"],
        "m_data_parameters": ["ireallut", "i1", "default_str_len"],
        "m_optprop_parameters": ["stddev_atol"],
    },
    "tests/test_geometry/test_geometry.F90": {
        "m_data_parameters": ["mpiint", "ireals", "one", "zero", "default_str_len"],
    },
    "tests/test_helper_functions/test_helper_functions.F90": {
        "m_helper_functions": [
            "compute_normal_3d", "cstr", "determine_normal_direction",
            "distance_to_edge", "imp_reduce_mean",
        ],
    },
    "tests/test_LUT_3_10/test_LUT_3_10.F90": {"m_data_parameters": ["i1"]},
    "tests/test_LUT_3_24/test_LUT_3_24.F90": {
        "m_boxmc": ["t_boxmc"],
        "m_data_parameters": ["i1"],
    },
    "tests/test_LUT_8_10/test_LUT_8_10.F90": {"m_data_parameters": ["i1"]},
    "tests/test_mmap/test_mmap.F90": {"m_mmap": ["binary_file_to_mmap"]},
    "tests/test_netcdfio/test_netcdfio.F90": {"m_data_parameters": ["iintegers"]},
    "tests/test_petsc_scatterToZero/test_petsc_scatterToZero.F90": {
        "m_data_parameters": ["zero", "one"],
    },
    "tests/test_petsc_sort/test_petsc_sort.F90": {
        "m_data_parameters": ["ireals", "zero", "one"],
        "m_helper_functions": ["chkerr"],
    },
    "tests/test_plexrt_fish/test_plexrt_fish.F90": {"m_helper_functions": ["chkerr"]},
    "tests/test_plexrt_nca/test_plexrt_nca.F90": {"m_data_parameters": ["default_str_len"]},
    "tests/test_plexrt_rayli/test_plexrt_rayli.F90": {"m_helper_functions": ["chkerr"]},
    "tests/test_plexrt_rrtmg_lw_sw/test_plexrt_rrtmg_lw_sw.F90": {
        "m_data_parameters": ["default_str_len"],
        "m_helper_functions": ["chkerr"],
    },
    "tests/test_pprts_absorption_by_coeff_divergence/test_pprts_absorption_by_coeff_divergence.F90": {
        "m_helper_functions": ["deg2rad", "imp_allreduce_mean"],
    },
    "tests/test_pprts_rrtm_icollapse/test_pprts_rrtm_icollapse.F90": {
        "iso_fortran_env": ["real32", "real64"],
        "m_data_parameters": ["one"],
    },
    "tests/test_pprts_slope_correction/test_pprts_slope_correction.F90": {
        "m_data_parameters": ["ireallut", "pi", "pi64"],
        "m_tenstream_options": ["read_commandline_options"],
        "m_dyn_atm_to_rrtmg": ["destroy_tenstr_atm", "print_tenstr_atm"],
    },
    "tests/test_pprts_solution_vecscale/test_pprts_solution_vecscale.F90": {
        "m_data_parameters": ["init_mpi_data_parameters", "ireallut", "zero", "pi"],
        "m_tenstream_options": ["read_commandline_options"],
    },
    "tests/test_pprts_specint/test_pprts_specint.F90": {
        "iso_fortran_env": ["real32", "real64"],
        "m_data_parameters": ["one", "share_dir"],
    },
    "tests/test_pprts_symmetry/test_pprts_symmetry.F90": {
        "m_data_parameters": ["pi"],
        "m_optprop": ["t_optprop"],
    },
    "tests/test_ranlux/test_ranlux.F90": {"m_data_parameters": ["ireals"]},
    "tests/test_rrtm_lw_Bsrfc/test_rrtm_lw_Bsrfc.F90": {
        "iso_fortran_env": ["real32", "real64"],
    },
    "tests/test_schwarzschild/test_schwarzschild.F90": {"m_schwarzschild": ["schwarzschild"]},
    "tests/test_search/test_search.F90": {
        "m_data_parameters": ["mpiint", "init_mpi_data_parameters"],
    },
    "tests/test_tenstr_atm/test_tenstr_atm.F90": {
        "iso_fortran_env": ["real32", "real64"],
    },
    "tests/test_twostr/test_twostr.F90": {"m_data_parameters": ["pi"]},
    "tests/test_wedge_boxmc/test_wedge_boxmc.F90": {
        "m_boxmc": ["t_boxmc"],
        "m_data_parameters": ["i1"],
        "m_optprop_parameters": ["stddev_atol"],
    },
    "tests/test_wedge_boxmc_18_8/test_wedge_boxmc_18_8.F90": {
        "m_boxmc": ["t_boxmc"],
        "m_data_parameters": ["i1"],
        "m_optprop_parameters": ["stddev_atol"],
    },
    "tests/test_wedge_boxmc_5_8/test_wedge_boxmc_5_8.F90": {
        "m_boxmc": ["t_boxmc"],
        "m_data_parameters": ["i1"],
        "m_optprop_parameters": ["stddev_atol"],
    },
    "tests/test_wedge_boxmc_5_8_spherical/test_wedge_boxmc_5_8_spherical.F90": {
        "m_boxmc": ["t_boxmc"],
        "m_data_parameters": ["i1"],
        "m_optprop_parameters": ["stddev_atol"],
        "m_helper_functions": ["itoa"],
        "m_boxmc_geometry": ["setup_default_unit_wedge_geometry"],
    },
    "tests/test_wedge_param_phi/test_wedge_param_phi.F90": {
        "m_boxmc": ["t_boxmc_wedge_18_8"],
        "m_data_parameters": ["ireals", "ireallut", "default_str_len", "i1", "i2", "i3", "i4", "i5"],
        "m_optprop": ["t_optprop_wedge_18_8"],
        "m_tenstream_options": ["read_commandline_options"],
        "m_helper_functions": ["rmse", "chkerr", "get_arg", "itoa", "ind_nd_to_1d", "ind_1d_to_nd"],
        "m_search": ["find_real_location"],
    },
}


def join_continuation_lines(lines):
    """Join Fortran continuation lines (ending with &) into logical lines.
    Returns list of (logical_line, [physical_line_indices]).
    """
    result = []
    i = 0
    while i < len(lines):
        line = lines[i]
        physical = [i]
        # Check if line continues (strip comments first for the check)
        stripped = line.rstrip()
        while stripped.endswith('&'):
            i += 1
            if i >= len(lines):
                break
            physical.append(i)
            stripped = lines[i].rstrip()
        result.append((''.join(lines[j] for j in physical), physical))
        i += 1
    return result


def remove_symbols_from_use(logical_line, module_name_lower, symbols_lower):
    """Remove specified symbols from a 'use module, only: ...' statement.

    Returns (new_logical_line_or_None, removed_count).
    new_logical_line_or_None is None if the entire use statement should be removed.
    """
    # Match: use <whitespace> <module> <whitespace> , <whitespace> only: <stuff>
    # Case-insensitive
    pat = re.compile(
        r'^(\s*use\s+)(' + re.escape(module_name_lower) + r')(\s*,\s*only\s*:\s*)(.+)$',
        re.IGNORECASE | re.DOTALL
    )
    m = pat.match(logical_line.strip('\n'))
    if not m:
        return logical_line, 0

    prefix = m.group(1)
    mod = m.group(2)
    only_part = m.group(3)
    sym_text = m.group(4)

    # Collect the indentation of the original first physical line
    indent = re.match(r'(\s*)', logical_line).group(1)

    # Parse symbol list: split on commas, handle & continuations
    sym_text_clean = re.sub(r'&\s*\n\s*', '', sym_text)
    sym_text_clean = sym_text_clean.rstrip()

    # Split on commas, stripping each symbol
    raw_syms = [s.strip() for s in sym_text_clean.split(',')]
    raw_syms = [s for s in raw_syms if s]  # drop empty

    # Remove the specified symbols (case-insensitive, but preserve original case)
    kept = [s for s in raw_syms if s.lower() not in symbols_lower]
    removed = len(raw_syms) - len(kept)

    if removed == 0:
        return logical_line, 0

    if not kept:
        # Remove entire use statement
        return None, removed

    # Rebuild as single line (fprettify will reformat)
    new_line = indent + 'use ' + mod + ', only: ' + ', '.join(kept) + '\n'
    return new_line, removed


def process_file(filepath, removals_for_file):
    """Apply removals to a file. Returns number of symbols removed."""
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"  SKIP (not found): {filepath}", file=sys.stderr)
        return 0

    total_removed = 0
    new_lines = list(lines)
    # Track which physical lines have been processed/deleted
    deleted = set()

    # Process logical lines to find use statements
    i = 0
    while i < len(new_lines):
        if i in deleted:
            i += 1
            continue

        line = new_lines[i]
        stripped = line.strip().lower()

        # Quick filter: must start with 'use '
        if not stripped.startswith('use '):
            i += 1
            continue

        # Collect continuation lines
        physical = [i]
        j = i
        while new_lines[j].rstrip().endswith('&'):
            j += 1
            if j >= len(new_lines):
                break
            physical.append(j)

        logical = ''.join(new_lines[k] for k in physical)

        # Determine which module this is
        use_pat = re.compile(r'^\s*use\s+(\w+)', re.IGNORECASE)
        um = use_pat.match(logical)
        if not um:
            i += 1
            continue

        mod_lower = um.group(1).lower()
        if mod_lower not in removals_for_file:
            i = physical[-1] + 1
            continue

        syms_to_remove = set(removals_for_file[mod_lower])
        new_logical, n = remove_symbols_from_use(logical, mod_lower, syms_to_remove)
        total_removed += n

        if n > 0:
            if new_logical is None:
                # Delete all physical lines for this use statement
                for k in physical:
                    deleted.add(k)
            else:
                # Replace first physical line, delete the rest
                new_lines[physical[0]] = new_logical
                for k in physical[1:]:
                    deleted.add(k)

        i = physical[-1] + 1

    if total_removed == 0:
        return 0

    # Write result
    with open(filepath, 'w') as f:
        for k, line in enumerate(new_lines):
            if k not in deleted:
                f.write(line)

    return total_removed


def main():
    root = Path(__file__).parent.parent  # tenstream root

    grand_total = 0
    for rel_path, mod_removals in sorted(REMOVALS.items()):
        filepath = root / rel_path
        n = process_file(str(filepath), mod_removals)
        if n > 0:
            print(f"  {rel_path}: removed {n} symbol(s)")
            grand_total += n
        elif not filepath.exists():
            pass  # already warned
        else:
            print(f"  {rel_path}: nothing removed (symbols may already be gone or not found)")

    print(f"\nTotal symbols removed: {grand_total}")


if __name__ == '__main__':
    main()
