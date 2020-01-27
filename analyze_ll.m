
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
};
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

list_methods = {
    '_1s2ca';
    '_1s2c';
    '_2s3ci';
    '_2s3ct';
};

n = size(list_species, 1);

m = size(list_methods, 1);

for i = 1:n
    species = list_species{i};

    species_folder = ['test_search/est_results/',species];
    param_folder = [species_folder,'/params/'];
    ll_folder = [species_folder,'/ll/'];
    if ~exist(ll_folder)
        mkdir(ll_folder)
    end
    mean_folder = [species_folder,'/mean/'];
    if ~exist(mean_folder)
        mkdir(mean_folder)
    end
    for j = 1:m
        method = list_methods{j};

        load(['test_search/matdata/',species,'_data.mat'])
        S = omat';
        s1 = S(1,:);
        s1 = s1(s1~=0);
        s2 = S(2,:);
        s2 = s2(s2~=0);
        
        run(['ll',method])
    end
end
