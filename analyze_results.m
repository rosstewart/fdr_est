
clear

% list_species = {
%     'H.sapiens3';
%     'M.musculus2';
%     'M.musculus3';
%     };
% clear

list_species = {
'A.thaliana'
% 'C.elegans'
'D.melanogaster'
'E.coli'
'H.sapiens2'
'H.sapiens3'
'M.musculus'
'M.musculus2'
'M.musculus3'
% 'S.cerevisiae'
'S.cerevisiae2'
'S.cerevisiae3'
}
% species = 'M.musculus'
% species = 'H.sapiens'
% species = 'H.sapiens2'
% species = 'H.sapiens3'
% species = 'H.sapiens4'
% species = 'C.elegans'
% species = 'D.melanogaster'
% species = 'S.cerevisiae'
% species = 'S.cerevisiae2'
% species = 'S.cerevisiae3'
% species = 'E.coli'
% species = 'A.thaliana'

% species

n = size(list_species, 1)
for i = 1:n
    species = list_species{i};
%     run parse_psm.m
    load(['test_search/matdata/',species,'_data.mat'])
    run_all
end
