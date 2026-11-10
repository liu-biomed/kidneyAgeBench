#!/bin/bash
# Run the full reanalysis for one cohort set: KSET=original|corrected bash run_all.sh
S="$(cd "$(dirname "$0")" && pwd)"; L="$S/../logs"; cd /tmp
for s in 02a_reproduce 02b_gtex 02c_loco 02d_sensitivity 02e_biology 02e2_arm_composition 02f_clinical 02g_drug; do
  echo "### $s" ; Rscript "$S/$s.R" > "$L/${s}_${KSET}.log" 2>&1 || { echo "FAILED $s"; tail -20 "$L/${s}_${KSET}.log"; exit 1; }
done; echo ALL_OK
