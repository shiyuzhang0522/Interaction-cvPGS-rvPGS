#!/bin/bash
# ------------------------------------------------------------------
# Calculate PRS with PLINK2 based on LDpred2 and LASSOsum2 estimates
# Author: Shelley
# Date: 2025-08-01
# Usage: bash 06_run_LDpred2_LASSOsum2_plink2.sh <CV_number> <chr>
# ------------------------------------------------------------------

set -euo pipefail

CV_number=$1   # e.g., 01
chr=$2         # e.g., 7

# ---- Software ----
PLINK2=plink2   # assumes plink2 is in PATH

# ---- Genotype ----
BEDFILE="chr${chr}.imputation.reset.varid"

# ---- Sample keep lists ----
for FILE in \
    "cv${CV_number}_train.txt" \
    "cv${CV_number}_tune.txt" \
    "cv${CV_number}_val.txt"; do

  KEEP_FILE="${FILE}.keep"
  if [ -f "$KEEP_FILE" ]; then
    echo "[NOTICE] $KEEP_FILE already exists, using existing file."
  else
    awk '{print $1, $1}' "$FILE" > "$KEEP_FILE"
    orig_lines=$(wc -l < "$FILE")
    keep_lines=$(wc -l < "$KEEP_FILE")
    if [ "$orig_lines" -ne "$keep_lines" ]; then
      echo "[ERROR] Line count mismatch for $FILE and $KEEP_FILE ($orig_lines vs $keep_lines)"
      exit 1
    else
      echo "[OK] $FILE and $KEEP_FILE have matching lines: $orig_lines"
    fi
  fi
done

TRAIN="cv${CV_number}_train.txt.keep"
TUNE="cv${CV_number}_tune.txt.keep"
VALI="cv${CV_number}_val.txt.keep"

# ---- Output directories ----
LD_OUT="PRS/LDpred2/${CV_number}/"
LS_OUT="PRS/LASSOsum2/${CV_number}/"
mkdir -p "$LD_OUT" "$LS_OUT"

# ---- LDpred2 PRS ----
LD_BETA="LDpred2_CV${CV_number}_prs_plink2.txt"

N=$(awk '{print NF; exit}' "$LD_BETA")
BETA_COL="7-${N}"

$PLINK2 --score "$LD_BETA" cols=+scoresums,-scoreavgs header no-mean-imputation \
    --score-col-nums $BETA_COL \
    --bfile "$BEDFILE" --keep "$TRAIN" --threads 4 \
    --out "${LD_OUT}/CV${CV_number}_chr${chr}_prs_train"

$PLINK2 --score "$LD_BETA" cols=+scoresums,-scoreavgs header no-mean-imputation \
    --score-col-nums $BETA_COL \
    --bfile "$BEDFILE" --keep "$TUNE" --threads 4 \
    --out "${LD_OUT}/CV${CV_number}_chr${chr}_prs_tune"

$PLINK2 --score "$LD_BETA" cols=+scoresums,-scoreavgs header no-mean-imputation \
    --score-col-nums $BETA_COL \
    --bfile "$BEDFILE" --keep "$VALI" --threads 4 \
    --out "${LD_OUT}/CV${CV_number}_chr${chr}_prs_validation"

# ---- LASSOsum2 PRS ----
LS_BETA="LASSOsum2_CV${CV_number}_prs_plink2.txt"

N=$(awk '{print NF; exit}' "$LS_BETA")
BETA_COL="7-${N}"

$PLINK2 --score "$LS_BETA" cols=+scoresums,-scoreavgs header no-mean-imputation \
    --score-col-nums $BETA_COL \
    --bfile "$BEDFILE" --keep "$TRAIN" --threads 4 \
    --out "${LS_OUT}/CV${CV_number}_chr${chr}_prs_train"

$PLINK2 --score "$LS_BETA" cols=+scoresums,-scoreavgs header no-mean-imputation \
    --score-col-nums $BETA_COL \
    --bfile "$BEDFILE" --keep "$TUNE" --threads 4 \
    --out "${LS_OUT}/CV${CV_number}_chr${chr}_prs_tune"

$PLINK2 --score "$LS_BETA" cols=+scoresums,-scoreavgs header no-mean-imputation \
    --score-col-nums $BETA_COL \
    --bfile "$BEDFILE" --keep "$VALI" --threads 4 \
    --out "${LS_OUT}/CV${CV_number}_chr${chr}_prs_validation"

echo "PRS calculation completed for CV${CV_number}, chr${chr}"
echo "Script finished at: $(date)"
