#!/bin/bash

#for psm in test_search/HeLa/HeLa*_nod.tsv;
for psm in test_search/nist/*_hcd_nod.tsv;
do
	psm=$(basename $psm)
	sp=${psm/_nod.tsv/}
	echo $sp
	#nohup matlab -nodisplay -r "parse_psm_fn('test_search/HeLa/', '$sp')" > log/parsing_${sp}.log 2>&1 &
	nohup matlab -nodisplay -r "parse_psm_fn('test_search/nist/', '$sp')" > log/parsing_${sp}.log 2>&1 &
done




