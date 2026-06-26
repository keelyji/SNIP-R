#!/usr/bin/env bash
#
# tests/run_demo.sh
#
# End-to-end smoke test for the SNIP-R pipeline. Runs each CLI script
# against the bundled example data and verifies the expected outputs
# appear. Intended for a person evaluating whether the pipeline works
# on their machine ("can I install this and reproduce a known result?").
#
# Usage:
#   bash tests/run_demo.sh
#
# Or inside the Docker image:
#   docker run --rm -v "$PWD":/work -w /work snipr:latest bash tests/run_demo.sh

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

OUT="${SNIPR_DEMO_OUT:-$ROOT/demo_output}"
rm -rf "$OUT"
mkdir -p "$OUT"
echo "[demo] output dir: $OUT"

# ---- Part B step 1 --------------------------------------------------------
echo
echo "[demo] === Part B: prep flanks ==="
Rscript scripts/partB_prep_flanks.R \
    --regions example_data/Finaldf.bed.txt \
    --out-dir "$OUT" \
    --skip-fasta

# We cannot run the CRISPick web step in CI. Use the bundled CRISPick
# outputs instead, which correspond to a previously-run CRISPick session
# on these flanks (Feb 2026).
echo
echo "[demo] === Part B: pair gRNAs (using bundled CRISPick outputs) ==="
Rscript scripts/partB_pair_grnas.R \
    --left-crispick example_data/CRISPick/Left-sgrna-designs.txt \
    --right-crispick example_data/CRISPick/Right-sgrna-designs.txt \
    --conversion-table "$OUT/ROItoFlank.tsv" \
    --out "$OUT/paired_grnas.tsv"

# Sanity-check
n_pairs=$(($(wc -l < "$OUT/paired_grnas.tsv") - 1))
n_rois=$(awk -F'\t' 'NR>1{print $1}' "$OUT/paired_grnas.tsv" | sort -u | wc -l)
echo "[demo] paired_grnas.tsv: $n_pairs pairs across $n_rois ROIs"
[ "$n_pairs" -ge 100 ] || { echo "[demo] FAIL: too few pairs"; exit 1; }

# ---- Part C ---------------------------------------------------------------
echo
echo "[demo] === Part C: screen analysis ==="
Rscript scripts/partC_screen_analysis.R \
    --counts example_data/SupplementaryTable4_RawScreenData.xlsx \
    --config config/example_partC_config.yaml \
    --out-dir "$OUT"

# Required outputs (pandonor_*.tsv use canonical condition names; volcano PDFs
# use the pandonor.display_names map, so volcano_Acute.pdf / volcano_Chronic.pdf).
required=(
    pandonor_Acute_Low.tsv
    pandonor_Chronic_Neg.tsv
    per_pair_zscores_donor1.tsv
    per_region_zscores_donor1.tsv
    volcano_Acute.pdf
    volcano_Chronic.pdf
    concordance_Acute_Low.pdf
)
for f in "${required[@]}"; do
    [ -s "$OUT/$f" ] || { echo "[demo] FAIL: missing or empty $f"; exit 1; }
done

# IFNG-locus sanity check: SimpleName "0" (IFNG promoter) should be a
# strong hit in the Acute_Low pan-donor analysis (mean Z > 1).
ifng_meanZ=$(awk -F'\t' '$1=="0" {print $3}' "$OUT/pandonor_Acute_Low.tsv")
echo "[demo] IFNG promoter (SimpleName=0) Acute_Low mean Z = $ifng_meanZ"
awk -v z="$ifng_meanZ" 'BEGIN{exit !(z+0 > 1)}' \
    || { echo "[demo] FAIL: IFNG promoter mean Z below expected threshold"; exit 1; }

echo
echo "[demo] === ALL CHECKS PASSED ==="
echo
echo "[demo] Artifacts written to: $OUT"
echo "[demo] Files produced:"
ls -1 "$OUT" | sed 's/^/[demo]   /'
echo
echo "[demo] Try opening one of the PDFs (e.g. volcano_Acute_Low.pdf) to see"
echo "[demo] the analysis output, or the pandonor_Acute_Low.tsv table to see"
echo "[demo] the per-region Z-scores and t-test p-values."
