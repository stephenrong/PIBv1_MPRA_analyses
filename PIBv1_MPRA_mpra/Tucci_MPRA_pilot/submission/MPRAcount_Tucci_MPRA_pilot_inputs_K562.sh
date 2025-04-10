#!/bin/sh
conda activate workflow-mpra-oligo-barcode
sbatch -t 24:00:00 --mem=255GB -c 12 --partition=bigmem --wrap "cromwell run /gpfs/ysm/project/reilly/sr2446/Projects/TucciMPRAPilot/workflows/MPRAoligo/MPRAcount.wdl --inputs MPRAcount_Tucci_MPRA_pilot_inputs_K562.json"
