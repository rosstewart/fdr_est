#!/bin/bash

# Create scores directory if it doesn't exist
mkdir -p ./scores

# Loop over all .mgf files in the ./mgf_splits/ directory
for mgf_file in ./mgf_splits_7-16_seq_length/*.mgf; do
    # Extract the filename without the .mgf extension
    filename=$(basename "$mgf_file" .mgf)
    
    # Run pepnovo on the current .mgf file and output results to a corresponding file in ./scores
    pepnovo -file "$mgf_file" -model CID_IT_TRYP -PTMs C+57:M+16 -digest TRYPSIN -correct_pm -num_solutions 10 > "./scores_trypsin_correct_pm_7-16/pepnovo_scores_${filename}.txt" &
done

