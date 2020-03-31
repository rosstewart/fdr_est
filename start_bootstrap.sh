#!/bin/bash

bootratio=1

for bi in `seq 1 20 200`;
do
	nohup matlab -nodisplay -r "bootstrap_all_fn(20, $bi, $bootratio)" > matlab${bi}.log 2>&1 &
done




