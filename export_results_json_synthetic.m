
clear

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

%json_dir = 'test_search/est_results_nist/json/';
json_dir = 'synthetic/fdr_result/json/';
if ~exist(json_dir)
    mkdir(json_dir)
end

% method = '_1s2c';

%results_folder = 'test_search/est_results_nist/';
results_folder = 'synthetic/fdr_result/';


list_methods = {
    '_1s2ca';
    '_1s2c';
    '_2s3ci';
%     '_2s3ct';
};

n_sp = size(list_species, 1);

n_m = size(list_methods, 1);

results_arr = {};

for sp_i = 1:n_sp
    species = list_species{sp_i}
    load(['synthetic/matdata/',species,'_data.mat'])
    species_folder = [results_folder,species];
    
    for me_j = 1:n_m
%         run parse_psm.m
        method = list_methods{me_j};
        method;

        load([species_folder, '/params/', method, '.mat'])

        if strcmp(method, '_1s2c')
            algo = '1S2D skew normal';
            alpha = theta.alpha;
            u_c = theta.theta_c.u;
            sigma_c = theta.theta_c.sigma;
            lambda_c = theta.theta_c.lambda;
            u_i = theta.theta_i.u;
            sigma_i = theta.theta_i.sigma;
            lambda_i = theta.theta_i.lambda;
        elseif strcmp(method, '_1s2ca')
            algo = '1S2D gamma & gaussian';
            alpha = theta.alpha;
            u_c = theta.theta_c.u;
            sigma_c = theta.theta_c.sigma;
            a_i = theta.theta_i.a;
            b_i = theta.theta_i.b;
            gamma_i = theta.theta_i.gamma;
        elseif strcmp(method, '_2s3ci') || strcmp(method, '_2s3ct')
            algo = '2S3D skew normal';
            if strcmp(method, '_2s3ct')
                algo = [algo, ' truncated'];
            end
            alpha = theta.alpha;
            beta = theta.beta;
            [u_c, sigma_c, lambda_c] = unpack_skntheta(theta.theta_c);
            [u_i, sigma_i, lambda_i] = unpack_skntheta(theta.theta_i);
            [u_i2, sigma_i2, lambda_i2] = unpack_skntheta(theta.theta_i2);
        end

        run(['result_curv', method])

        % y1
        % yi1
        % yc

        obj = {};
        obj.ds = species;
        obj.algo = algo;
        obj.s1 = s1;
        obj.pdfM = y1;
        obj.pdfI1 = yi1;
        obj.pdfC = yc;
        obj.cdfTrue = h1emp;
        obj.cdf = h1;
        obj.fdr = fdr_curv;
        obj.t1p = thres;
        obj.t01p = thres01;
        obj.t10p = thres10;
        obj.ll = ll_1;
        obj.deltaCdf = sdcdf;
        obj.pi_C = alpha;
        
        % obj.t1p
        figure
        hold on;
        histogram(s1, 'Normalization', 'pdf')
        plot(s1, y1)
        plot(s1, yi1)
        plot(s1, yc)
        plot(s1, h1emp)
        plot(s1, h1)
        plot(s1, fdr_curv)
        results_arr{end+1} = obj;
    end
    json = jsonencode(results_arr);
%     json
    fid = fopen([json_dir, species, '.json'],'wt');
    fprintf(fid, json);
    fclose(fid);
    results_arr = {};
end


