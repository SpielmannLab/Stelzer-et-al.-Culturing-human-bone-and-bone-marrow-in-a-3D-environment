#! /bin/bash
# Run this script to carry out the entire Single-Cell-RNA analysis with spliced and unspliced data added
# STEP1. Make sure you have run the appropriate pipelines in the SpielmannLab/scPreprocess repository
# STEP2. Provide all the details in the file "sc_p_sample_params.yaml"
# STEP3. Close this file and execute it using ./sc_p_sample_main.sh

sbatch sc_p_sample_core.sh
