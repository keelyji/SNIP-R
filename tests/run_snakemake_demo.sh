#!/usr/bin/env bash
#
# tests/run_snakemake_demo.sh
#
# End-to-end smoke test for the SNIP-R Snakemake workflow. Exercises all
# three rules (prep_flanks, pair_grnas, screen_analysis) end-to-end against
# the bundled example data.
#
# Uses skip_fasta: true to bypass the bedtoolsr/bedtools FASTA-extraction
# step in prep_flanks; this matches what tests/run_demo.sh does, and means
# the test runs without bedtools / bedtoolsr being installed. The goal here
# is to verify the Snakemake workflow itself (rule chaining, config wiring,
# file dependencies), not to re-test the generic bedtools wrapper.
#
# If you want to ALSO exercise the FASTA-extraction step, decompress
# notebooks/PartB_SNIP-R_sgRNApairing/chr12.fa.gz to example_data/chr12.fa,
# add `reference: example_data/chr12.fa` to the config below, and remove the
# `skip_fasta: true` line. bedtoolsr (and a bedtools binary) must be on PATH.
# This is already the case inside the SNIP-R Docker image.
#
# Usage:
#   bash tests/run_snakemake_demo.sh
#
# Or inside the Docker image:
#   docker run --rm -v "$PWD":/work -w /work snipr:latest \
#       bash tests/run_snakemake_demo.sh

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

# ---- stage / cleanup ------------------------------------------------------
STAGED_CONFIG="config.yaml"
OUT="snakemake_output"

# Cleanup behaviour depends on whether the test succeeded:
#  - On SUCCESS: keep $OUT so the user can inspect outputs (volcanos, tables);
#                only sweep the test-generated config.yaml and .snakemake/.
#  - On FAILURE: remove $OUT too so the next run starts from a clean slate.
SUCCESS=0
cleanup() {
    if [ "$SUCCESS" = "1" ]; then
        rm -f "$STAGED_CONFIG"
        rm -rf .snakemake
    else
        rm -f "$STAGED_CONFIG"
        rm -rf "$OUT" .snakemake
    fi
}
trap cleanup EXIT

echo "[snakemake-demo] writing $STAGED_CONFIG"
cat > "$STAGED_CONFIG" <<EOF
regions:        example_data/Finaldf.bed.txt
crispick_left:  example_data/CRISPick/Left-sgrna-designs.txt
crispick_right: example_data/CRISPick/Right-sgrna-designs.txt
counts:         example_data/SupplementaryTable4_RawScreenData.xlsx
screen_config:  config/example_partC_config.yaml
out_dir:        $OUT
skip_fasta:     true
EOF

# ---- workflow -------------------------------------------------------------
echo
echo "[snakemake-demo] === dry run ==="
snakemake -n -j1 all

echo
echo "[snakemake-demo] === full run ==="
snakemake -j1 all

# ---- verify ---------------------------------------------------------------
echo
echo "[snakemake-demo] === verifying outputs ==="

required=(
    "$OUT/ROItoFlank.tsv"
    "$OUT/paired_grnas.tsv"
    "$OUT/pandonor_Acute_Low.tsv"
    "$OUT/pandonor_Chronic_Neg.tsv"
    "$OUT/volcano_Acute.pdf"
    "$OUT/volcano_Chronic.pdf"
    # partC_done.flag is intentionally a 0-byte sentinel created by
    # the screen_analysis rule; Snakemakes "100% done" report is its own
    # verification. The data files below are the real check.
)
for f in "${required[@]}"; do
    [ -s "$f" ] || { echo "[snakemake-demo] FAIL: missing or empty $f"; exit 1; }
done

# Same IFNG-promoter sanity check used by tests/run_demo.sh
ifng_meanZ=$(awk -F'\t' '$1=="0" {print $3}' "$OUT/pandonor_Acute_Low.tsv")
echo "[snakemake-demo] IFNG promoter (SimpleName=0) Acute_Low mean Z = $ifng_meanZ"
awk -v z="$ifng_meanZ" 'BEGIN{exit !(z+0 > 1)}' \
    || { echo "[snakemake-demo] FAIL: IFNG promoter mean Z below expected threshold"; exit 1; }

echo
echo "[snakemake-demo] === ALL CHECKS PASSED ==="
echo "[snakemake-demo] Output files (preserved for inspection) in: $OUT/"
echo "[snakemake-demo]   try opening $OUT/volcano_Chronic.pdf or $OUT/pandonor_Acute_Low.tsv"
echo "[snakemake-demo] Temp artifacts ($STAGED_CONFIG, .snakemake/) will be removed on exit."
SUCCESS=1
