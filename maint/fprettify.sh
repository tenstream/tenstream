#!/bin/bash
set -eu -o pipefail

function show_help {
echo "script to reformat the codebase with fprettify"
echo "  Usage: misc/fprettify.sh [OPTION]"
echo ""
echo "  -x=*, --fprettify={path}         Set path to fprettify. If not present, we try to download it"
echo "  -n, --dry-run                    Only check for illegal lines, dont change anything"
echo "  -d, --diff                       print diff to stdout, dont change anything"
echo "  -v, --verbose                    Increase verbosity"
echo "  -f, --files                      Specify files on which to run"
echo "  -h, --help                       Show this help text and exit"
echo "  -p, --parallel <ncpu>            set the number of parallel threads or we use as many as we find in /proc/cpuinfo "
echo "  -s, --staging                    only look at files that are currently changed "
}

SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
PROJECT_ROOT=$(readlink -f "$SCRIPTDIR/../")
ERRLOG=".fprettify.errlog"

cd $PROJECT_ROOT
rm -f $ERRLOG
touch $ERRLOG

PYTHON=$(which "python3")
BIN=$PROJECT_ROOT/external/fprettify/fprettify.py
if [[ ! -e $BIN ]]; then
  git clone https://github.com/pseewald/fprettify.git $(dirname $BIN)
  sed -i 's#COMMENT_LINE_STR = r"^!"#COMMENT_LINE_STR = r"^!|^@"#' $(dirname $BIN)/fprettify/fparse_utils.py
fi
DIFF=false
DRYRUN=false
VERBOSE=false
NAMEONLY=false

# Determine number of processors
if [ -e /proc/cpuinfo ]; then
  NCPU=$(grep -c ^processor /proc/cpuinfo)
else
  NCPU=4
fi

FILES_VIEW='default'
FILES="\
  src/*.fypp \
  src/*.F90 \
  src/*.inc \
  rrtmg/rrtmg/*.F90 \
  repwvl/*.F90 \
  specint/*.F90 \
  plexrt/*.F90 \
  c_wrapper/*.F90 \
  tests/**/*.F90 \
  examples/**/*.F90 \
  "


while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -p=*|--python=*)
    PYTHON="${key#*=}"
    shift # past argument=value
    ;;
    -x=*|--fprettify=*)
    BIN="${key#*=}"
    shift # past argument=value
    ;;
    -n|--dry-run)
    DRYRUN=true
    shift # past argument
    ;;
    -d|--diff)
    DIFF=true
    shift # past argument
    ;;
    -j|--parallel)
    NCPU="$2"
    shift # past argument
    shift # past value
    ;;
    -s|--staging)
    NAMEONLY=true
    FILES_VIEW="staging"
    shift # past argument
    ;;
    -f|--file)
    FILES="$2"
    FILES_VIEW="$2"
    shift # past argument
    shift # past value
    ;;
    -v|--verbose)
    VERBOSE=true
    shift # past argument
    ;;
    -h|--help)
    shift # past argument
    show_help
    exit 0
    ;;
    *)    # unknown option
			echo "unknown argument $key"
      show_help
      exit 1
    ;;
esac
done

# Print Arguments
if $VERBOSE ; then
	echo "Found arguments:"
	echo "--fprettify=$BIN"
	echo "--dry-run=$DRYRUN"
	echo "--diff=$DIFF"
	echo "--verbose=$VERBOSE"
	echo "--parallel=$NCPU"
	echo "--staging=$NAMEONLY"
  echo "--files=$FILES_VIEW"
	echo ""
fi

if $NAMEONLY; then FILES=$(git diff --name-only | egrep "\.inc$|\.F90$|\.fypp$"); fi

if [[ ! -e $PYTHON ]]; then
  echo "Could not find python3! ($PYTHON)"
  exit 1
fi

if [[ ! -e $BIN ]]; then
  echo "Could not find fprettify at $BIN"
  exit 2
fi

OPT="\
  --indent 2 \
  --line-length 132 \
  --whitespace 3 \
  --strict-indent \
  --enable-decl \
  --case 1 1 1 1 \
  "

if $DRYRUN; then OPT="$OPT --stdout"; fi
if $DIFF  ; then OPT="$OPT --diff"; fi

count=0

for f in $FILES
do
  (
  if $VERBOSE ; then echo "prettify => $f"; fi
  $PYTHON $BIN $OPT $f &> >(tee -a $ERRLOG)
  ) &
  count=$(( $count + 1 ))
  if [ $count -eq $NCPU ]; then
    count=0
    wait
  fi
done
wait
if $VERBOSE; then echo done; fi
exit $(wc -l < $ERRLOG)
