#!/bin/sh
set -eu
# Run from the example folder. Refuse to overwrite earlier results.
if [ -e outputs ]; then
    echo 'outputs already exists. Move it to another folder before repeating.' >&2
    exit 1
fi
mkdir outputs
./bin/rsbagel fc -i inputs/counts.tsv -o outputs/results -c 1
./bin/rsbagel bf -i outputs/results.foldchange -o outputs/results.bf -e inputs/essentials.tsv -n inputs/nonessentials.tsv -c 1,2 --seed 42 --workers 1
./bin/rsbagel pr -i outputs/results.bf -o outputs/results.pr -e inputs/essentials.tsv -n inputs/nonessentials.tsv
for name in results.foldchange results.normed_readcount results.bf results.pr; do
    cmp "expected/$name" "outputs/$name"
done
echo 'SUCCESS: all four scientific result files match expected outputs.'
echo 'Results are in outputs/results.pr and outputs/results.bf.'
