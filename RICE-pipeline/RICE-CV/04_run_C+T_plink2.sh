#!/bin/bash
# ------------------------------------------------------------------
# Calculate C+T based PRS using PLINK2 for train, tune, and validation sets
# Author: Shelley
# Date: 2025-08-03
# Usage: bash 04_run_C+T_plink2.sh <CV_number> <chr>
# Notice: please remember to update the variant id for the UKBB imputed data as chr:pos:ref:alt, not rsid!!!
# ------------------------------------------------------------------

set -euo pipefail

CV_number=$1   # e.g., 01
chr=$2         # e.g., 7

# ---- Software ----
PLINK2=plink2   # assumes plink2 is in PATH

# ---- Genotype ----
BEDFILE="chr${chr}.imputation.reset.varid"

# ---- Sample keep lists (FID IID, no header) ----
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

# ---- Output directory ----
C_T_OUT="C+T_PRS/${CV_number}/"
mkdir -p "$C_T_OUT"

# ---- PRS input files ----
Q_RANGE="${C_T_OUT}/C+T_q_range.txt"
P_VALUE="${C_T_OUT}/C+T_p_value.txt"
COEFF="${C_T_OUT}/C+T_prs_coeff.txt"

# ---- Run PLINK2 C+T PRS ----
echo "[INFO] Calculating C+T PRS for train set..."
$PLINK2 --q-score-range "$Q_RANGE" "$P_VALUE" header \
    --threads 2 \
    --score "$COEFF" cols=+scoresums,-scoreavgs header no-mean-imputation \
    --bfile "$BEDFILE" \
    --keep "$TRAIN" \
    --out "${C_T_OUT}/CV${CV_number}.chr${chr}.C+T_prs_train"

echo "[INFO] Calculating C+T PRS for tune set..."
$PLINK2 --q-score-range "$Q_RANGE" "$P_VALUE" header \
    --threads 2 \
    --score "$COEFF" cols=+scoresums,-scoreavgs header no-mean-imputation \
    --bfile "$BEDFILE" \
    --keep "$TUNE" \
    --out "${C_T_OUT}/CV${CV_number}.chr${chr}.C+T_prs_tune"

echo "[INFO] Calculating C+T PRS for validation set..."
$PLINK2 --q-score-range "$Q_RANGE" "$P_VALUE" header \
    --threads 2 \
    --score "$COEFF" cols=+scoresums,-scoreavgs header no-mean-imputation \
    --bfile "$BEDFILE" \
    --keep "$VALI" \
    --out "${C_T_OUT}/CV${CV_number}.chr${chr}.C+T_prs_validation"

echo "C+T PRS calculation completed for CV${CV_number}, chr${chr} at $(date)"
