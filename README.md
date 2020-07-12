# Decoy-Free FDR Estimation in Mass-spectrometry with Skew Normal Mixture Model

**Note:** Commands in documentation are assuming Linux environment. For other OS it will be similar but might need small modifications.

## Steps
1. Get peak data files (.mgf, .mzMl) using software or directly download from some database.
2. Run MSGF+ to search for matching peptides.
3. Use MSGF+ tools to convert .mzid result files into .tsv files.
4. Run parsing program to parse .tsv data and save into .mat format for MATLAB program.
5. Run program in MATLAB to estimate parameters and get FDR estimation.

## MSGF+ Searching Configuration
* Set decoy to off.
* Set number of matches to report as 10.
  * Only top 2 matches are used for methods reported in the paper, but considering further use top 10 to be reported is suggested.
* Matching parameters just set as normal.

## MSGF+ Command to Convert .mzid Results into .tsv Files
```sh
java -Xmx8G -cp msgf/MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv -showDecoy 1 -i <mzid_FILE>
```
change `<mzid_FILE>` to the path to the result file

## Searching Result Parsing
In order to run the program, results are parsed with python program and saved into .mat files for better performance.
1. Put the result files into a directory. Rename the files as `<dataset_id>_nod.tsv`. '_nod' stands for searching results without decoy.
2. Add config into the python code 'parse_msgf_psm.py'. Examples can be seen at top.
    * `data_source`: Switch to change between several set of results. This can be omitted if there is no such need.
    * `psm_dir`: Path of the directory holding .tsv files to be parse
    * `data_dir`: Place to hold the converted `.mat` files
3. After `parse_msgf_psm.py` is configured, run it with python.
```sh
python parse_msgf_psm.py
```

## Run FDR Estimation Program
`analyze_results.m` is the entrance script.
1. Configure the script:
    * `list_species`: A list of searching results to be analyzed. Each searching result is represented with `dataset_id`.
    Just put all `dataset_id` into the list. The word 'species' in the name is deprecated. 'datasets' would be more appropriate.
    * `data_dir`: Directory holding the saved `.mat` files resulted from previous parsing step.
    * `results_folder`: Place to hold results and plottings going to be create by the program.
2. Edit `run_all.m` file to choose the methods you want to run. (Please check our paper for details on each method.)
    * `run_1s2ca`: Gammar-gaussian model
    * `run_1s2c`: 1SMix, One sample skew normal mixture model.
    * `run_2s3ci`: 2SMix, Two samples skew normal mixture model.
3. Run the program in MATLAB:
    ```matlab
    >> analyze_results
    ```
