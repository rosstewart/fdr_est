#!/bin/bash

# Create scores directory if it doesn't exist
mkdir -p ./scores

# Loop over all .mgf files in the ./mgf_splits/ directory
for mgf_file in ./mgf_splits/*.mgf; do
    # Extract the filename without the .mgf extension
    filename=$(basename "$mgf_file" .mgf)
    
    # Run pepnovo on the current .mgf file and output results to a corresponding file in ./scores
    pepnovo -file "$mgf_file" -model CID_IT_TRYP -PTMs C+57:M+16 -digest NON_SPECIFIC -num_solutions 10 > "./scores/pepnovo_scores_${filename}.txt" &
done

