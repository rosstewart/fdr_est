function analyze_results_fn(species)
	species

	data_dir = 'test_search/matdata_nist/';
	results_folder = 'test_search/est_results_nist_allinits/';

    load(['test_search/matdata_nist/',species,'_data.mat'])
    run_all
end
