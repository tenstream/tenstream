#!/bin/bash
# Wrapper around build_dependencies.sh for CI jobs.
#
# PETSc + NetCDF are the slowest part of every CI job. The GitLab cache restores
# a previous build into $PETSC_DIR; this script reuses it if it was built with the
# same configuration (branch, precision, index size, compilers, MPI, build scripts)
# and is not older than PETSC_CACHE_MAX_AGE_DAYS. Otherwise it rebuilds from
# scratch and prunes the build intermediates so that the cache stays small.
#
# Env:
#   PETSC_CACHE_MAX_AGE_DAYS  rebuild if the cached build is older (default: 7)
#   PETSC_CACHE_REBUILD       set to non-empty to force a rebuild
set -euo pipefail

SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

: "${PETSC_DIR:?PETSC_DIR has to be set}"
: "${PETSC_ARCH:?PETSC_ARCH has to be set}"
MAX_AGE_DAYS=${PETSC_CACHE_MAX_AGE_DAYS:-7}
STAMP="$PETSC_DIR/.ci_build_stamp"

first_line() { "$@" 2>&1 | head -n 1 || true; }

signature() {
  echo "branch=${PETSC_BRANCH:-main} arch=$PETSC_ARCH precision=${PETSC_PRECISION:-} debugging=${PETSC_DEBUGGING:-} int64=${PETSC_64_INTEGERS:-}"
  echo "opts=${PETSC_OPTS:-}"
  env | grep -E '^PETSC_[A-Z]*FLAGS=' | sort || true
  echo "CC=${CC}: $(first_line ${CC} --version)"
  echo "CXX=${CXX}: $(first_line ${CXX} --version)"
  echo "FC=${FC}: $(first_line ${FC} --version)"
  echo "mpirun: $(first_line mpirun --version)"
  echo "scripts: $(cat "$SCRIPTDIR/build_dependencies.sh" "${BASH_SOURCE[0]}" | sha1sum | cut -d ' ' -f 1)"
}

SIG="$(signature)"

if [[ -f "$STAMP" && -z "${PETSC_CACHE_REBUILD:-}" ]]; then
  BUILT=$(head -n 1 "$STAMP")
  AGE_DAYS=$(( ($(date +%s) - BUILT) / 86400 ))
  if [[ "$(tail -n +2 "$STAMP")" != "$SIG" ]]; then
    echo "Cached PETSc build in $PETSC_DIR has a different configuration -- rebuilding:"
    diff <(tail -n +2 "$STAMP") <(echo "$SIG") || true
  elif (( AGE_DAYS >= MAX_AGE_DAYS )); then
    echo "Cached PETSc build in $PETSC_DIR is $AGE_DAYS days old (max $MAX_AGE_DAYS) -- rebuilding"
  elif [[ ! -e "$PETSC_DIR/$PETSC_ARCH/lib/pkgconfig/PETSc.pc" || ! -e "$PETSC_DIR/$PETSC_ARCH/lib/pkgconfig/netcdf-fortran.pc" ]]; then
    echo "Cached PETSc build in $PETSC_DIR is incomplete -- rebuilding"
  else
    echo "Reusing cached PETSc/NetCDF build in $PETSC_DIR (built $AGE_DAYS days ago)"
    echo "$SIG"
    exit 0
  fi
fi

# Only wipe/prune in GitLab CI: when run locally (e.g. misc/ci-docker-run.py),
# $PETSC_DIR may be a developer's checkout that we must not throw away.
if [[ "${GITLAB_CI:-}" == "true" ]]; then
  rm -rf "$PETSC_DIR"
  export PETSC_CLONE_OPTS=${PETSC_CLONE_OPTS:---depth 1}
fi

"$SCRIPTDIR/build_dependencies.sh"

if [[ "${GITLAB_CI:-}" == "true" ]]; then
  echo "Pruning PETSc build intermediates to keep the CI cache small"
  rm -rf "$PETSC_DIR/.git" "$PETSC_DIR/$PETSC_ARCH/obj" "$PETSC_DIR/$PETSC_ARCH/externalpackages"
  du -sh "$PETSC_DIR" || true
fi

{ date +%s; echo "$SIG"; } > "$STAMP"
