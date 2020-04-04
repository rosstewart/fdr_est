
clear

% list_species = {
%     'H.sapiens3';
%     'M.musculus2';
%     'M.musculus3';
%     };
% clear

list_species = {
% 'A.thaliana'
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
}
list_species = {
%     'c_elegans'
%     'drosophila'
%     'e_coli'
%     'human'
%     'mouse'
    'yeast'
}

data_dir = 'test_search/matdata_nist/';

results_folder = 'test_search/est_results_nist';

n = size(list_species, 1)
for i = 1:n
    species = list_species{i};
%     run parse_psm.m
    load([data_dir,species,'_data.mat'])
    run_all
end
