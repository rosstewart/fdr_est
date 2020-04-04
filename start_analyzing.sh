#!/bin/bash

for matfile in test_search/matdata_nist/*.mat
do
	matfile=$(basename $matfile)
	sp=${matfile/_data.mat/}
	echo $sp
	nohup matlab -nodisplay -r "analyze_results_fn('$sp')" > log/analyzing_${sp}.log 2>&1 &
done




