
clear

% list_species = {
%     'H.sapiens3';
%     'M.musculus2';
%     'M.musculus3';
%     };
% clear

%list_species = {
% 'a_thaliana'
% % 'C.elegans'
% 'D.melanogaster'
% 'E.coli'
% 'H.sapiens2'
% 'H.sapiens3'
% 'M.musculus'
% 'M.musculus2'
% 'M.musculus3'
% % 'S.cerevisiae'
% 'S.cerevisiae2'
% 'S.cerevisiae3'
%}
data_dir = 'synthetic/matdata/';

results_folder = 'synthetic/fdr_result/';
xlim_high = 65;
plotcdf = true;
plotthres = true;

% list_species = {
% 'HeLa01ng'
% 'HeLa1ng'
% 'HeLa10ng'
% 'HeLa50ng'
% 'HeLa100ng'
% 'HeLa01ng.2'
% 'HeLa1ng.2'
% 'HeLa10ng.2'
% 'HeLa50ng.2'
% 'HeLa100ng.2'
% 'HeLa01ng.3'
% 'HeLa1ng.3'
% 'HeLa10ng.3'
% 'HeLa50ng.3'
% 'HeLa100ng.3'
% }
% data_dir = 'test_search/matdata_hela/';
% 
% results_folder = 'test_search/est_results_hela';


% list_species = {
%     'c_elegans',
%     'h_sapiens',
%     'm_musculus',
%     's_cerevisiae'
% }

list_species = {
    'synthetic_1',
    'synthetic_2',
    'synthetic_3',
    'synthetic_4',
    'synthetic_5',
    'synthetic_6',
    'synthetic_7',
    'synthetic_8',
    'synthetic_9',
    'synthetic_10',
    'synthetic_11',
    'synthetic_12',
    'synthetic_13',
    'synthetic_14',
    'synthetic_15',
    'synthetic_16',
    'synthetic_17',
    'synthetic_18',
    'synthetic_19',
    'synthetic_20',
    'synthetic_21',
    'synthetic_22',
    'synthetic_23',
    'synthetic_24',
    'synthetic_25',
    'synthetic_26',
    'synthetic_27',
    'synthetic_28',
    'synthetic_29',
    'synthetic_30'
}
% 
% data_dir = 'test_search/matdata_nist/';
% 
% results_folder = 'test_search/est_results_nist';

figw = 240;
figh = 120;

c_color = [0.8500 0.3250 0.0980];
i1_color = [0.9290 0.6940 0.1250];
i2_color = [0 0.4470 0.7410];
mix_color = [0.4940 0.1840 0.5560];

cdfcolor = [0.6350 0.0780 0.1840];
legending = false;

legending2 = false;

linewidth = 1;

n = size(list_species, 1)
for i = 1:n
    % min_dcdf = 1;
    species = list_species{i};
%     run parse_psm.m
    load([data_dir,species,'_data' ...
        '.mat'])
    run_all
end
