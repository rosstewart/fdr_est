list_species = {
% 'M.musculus'
% 'M.musculus2'
% 'M.musculus3'
% 'H.sapiens2'
'H.sapiens3'
'C.elegans'
'D.melanogaster'
'S.cerevisiae'
'S.cerevisiae2'
'S.cerevisiae3'
'E.coli'
'A.thaliana'
}

n = size(list_species, 1)
for i = 1:n
    species = list_species{i};
    run parse_psm.m
end
