#!/bin/bash

for matfile in `ls test_search/matdata_hela/*.c*.mat`
do
	matfile=$(basename $matfile)
	sp=${matfile/_data.mat/}
	echo $sp
	nohup matlab -nodisplay -r "analyze_results_fn('$sp')" > log/analyzing_${sp}.log 2>&1 &
done




