#!/bin/bash

TMP=`sed -n ${SLURM_ARRAY_TASK_ID}p ../some_documents/158_gse_for_slurm_script.csv`
IFS='    ' read -r gse gpl <<< $TMP

Rscript start_diff_expression.R $gse $gpl
