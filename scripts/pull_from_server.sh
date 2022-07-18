#!/usr/bin/env bash
remote=${PROJECT_REMOTE_HOST:-TODO}
project=${REMOTE_PROJECT_DIR:-TODO}

echo "Pulling files from $remote:$project"

for pattern in $@
do
    rsync -hravz \
        --progress --partial --prune-empty-dirs --copy-dirlinks --relative \
        $remote:"$project/./$pattern" .
done
