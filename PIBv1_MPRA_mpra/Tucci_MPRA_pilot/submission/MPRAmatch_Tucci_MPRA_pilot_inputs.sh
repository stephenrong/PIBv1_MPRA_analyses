#!/bin/sh
conda activate workflow-mpra-oligo-barcode
sbatch -t 24:00:00 --mem=128GB -c 8 --wrap "cromwell run /gpfs/ysm/project/reilly/sr2446/Projects/TucciMPRAPilot/workflows/MPRAoligo/MPRAmatch.wdl --inputs MPRAmatch_Tucci_MPRA_pilot_inputs.json"
