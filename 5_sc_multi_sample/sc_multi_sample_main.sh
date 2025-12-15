#! /bin/bash
# Run this script to carry out the merging/integration of the multiple sample ending with de_gene identification
# STEP1. Make sure you have run the sc_p_sample_main.sh script
# STEP2. Provide all the details in the file "sc_multi_sample.yaml"
# STEP3. Close this file and execute it using ./sc_multi_sample_main.sh

sbatch sc_multi_sample_core.sh 
