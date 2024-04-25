#!/bin/bash
set -euo pipefail

PGDATA=/scratch-local/Fabian.Jakub/bmc_pgdata
mkdir -p $PGDATA

podman rm bmcdb || echo "no past db present"

podman run \
  --name bmcdb \
  -e POSTGRES_DB=bmcdb \
  -e POSTGRES_USER=bmc \
  -e POSTGRES_PASSWORD=bmcdb \
  -e PGDATA=/var/lib/postgresql/data/pgdata \
  -p 5432:5432 \
  -v $PGDATA:/var/lib/postgresql/data:z \
  postgres:16
