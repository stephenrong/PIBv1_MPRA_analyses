#!/bin/sh
conda activate workflow-mpra-oligo-barcode
sbatch -t 48:00:00 --mem=128GB -c 18 --partition=bigmem --wrap "cromwell run /gpfs/ysm/project/reilly/sr2446/Projects/TucciMPRAPilot/workflows/MPRAoligo/MPRAcount.wdl --inputs MPRAcount_Tucci_MPRA_pilot_inputs.json"
