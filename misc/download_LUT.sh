#!/bin/bash
set -eu -o pipefail

SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
PROJECT_ROOT=$(readlink -f $SCRIPTDIR/../)
DST=$PROJECT_ROOT/misc/LUT

MASK="${1:-LUT}"
DST="${2:-${DST}}"

echo "Looking for Tables with mask: $MASK and saving files to $DST"

mkdir -p $DST
cd $DST

CURL_BIN=$(command -v curl || true)
if [ -z $CURL_BIN ]; then
  echo "Could not find curl but I need it to download files"
  echo "Make sure we have curl available"
  exit 1
fi

BASEURL="https://www.meteo.physik.uni-muenchen.de/~Fabian.Jakub/TenstreamLUT/"
INDEX=$(curl -s $BASEURL)

RED="\033[0;31m"
GREEN="\033[0;32m"
NC='\033[0m' # No Color

for line in $INDEX
do
  if [[ ! -z $(echo $line | grep 'href=".*nc' ) ]]
  then
    fname=$(echo $line | cut -d'"' -f 2)
    if [[ "$fname" == *"$MASK"* ]]; then
      if [[ -e $DST/$fname ]]; then # file already exists
        # curl has a bug with the -C option if file already fully exists
        # as a work around: first check for file size and if it is already there, skip it
        REMOTE_FILE_LENGTH=$($CURL_BIN -I $BASEURL/$fname 2> /dev/null | grep Content-Length | awk '{print int($2)}')
        LOCAL_FILE_LENGTH=$(stat --printf="%s" $DST/$fname)
        if [ $REMOTE_FILE_LENGTH -eq $LOCAL_FILE_LENGTH ]; then
          printf "${GREEN}** Skipping Download of $fname because remote and local file sizes match: [${REMOTE_FILE_LENGTH} // ${LOCAL_FILE_LENGTH}]\n"
        else
          printf "${GREEN}** Resuming Download $fname${NC} [${REMOTE_FILE_LENGTH} // ${LOCAL_FILE_LENGTH}]\n"
            $CURL_BIN -C - --progress-bar $BASEURL/$fname --output $DST/$fname
        fi
      else
        printf "${GREEN}** Downloading $fname${NC}\n"
          $CURL_BIN --progress-bar $BASEURL/$fname --output $DST/$fname
      fi
    else
      printf "${RED}Not matching mask: $MASK -> Skipping $fname${NC}\n"
    fi
  fi
done

echo ""
echo "For convenience you may add the location of the luts to your global petscrc, i.e. run:"
echo ""
echo "  echo -lut_basename $DST/LUT >> $HOME/.petscrc"
echo ""
