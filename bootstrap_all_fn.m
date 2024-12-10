function bootstrap_all_fn(bootn, bootstart, bootratio)


list_species = {
%    'synthetic_1',
%    'synthetic_2',
%    'synthetic_3',
%    'synthetic_4',
%    'synthetic_5',
%    'synthetic_6',
%    'synthetic_7',
%    'synthetic_8',
%    'synthetic_9',
%    'synthetic_10',
%    'synthetic_11',
%    'synthetic_12',
%    'synthetic_13',
    'synthetic_14'
%    'synthetic_15',
%    'synthetic_16',
%    'synthetic_17',
%    'synthetic_18',
%    'synthetic_19',
%    'synthetic_20',
%    'synthetic_21',
%    'synthetic_22',
%    'synthetic_23',
%    'synthetic_24',
%    'synthetic_25',
%    'synthetic_26',
%    'synthetic_27',
%    'synthetic_28',
%    'synthetic_29',
%    'synthetic_30'
}


%list_species = {
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
% 
% 'HeLa01ng_2'
% 'HeLa1ng'
% 'HeLa10ng'
% 'HeLa50ng'
% 'HeLa100ng'

%'c_elegans'
%'h_sapiens'
%'m_musculus'
%'s_cerevisiae'

%'drosophila'
%'e_coli'
%'human'
%'mouse'

%}


n = size(list_species, 1)
for i = 1:n
    species = list_species{i}
%   run parse_synthetic_psm.py
    %load(['pepnovo/nist/matdata_7-16/',species,'_data.mat'])
    load(['synthetic/matdata/',species,'_data.mat'])
    bootstraping_fn(species, bootstart, bootn, bootratio)
end

end


