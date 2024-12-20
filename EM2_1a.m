function [alpha, u_c, sigma_c, a_i1, b_i1, gamma_i1] = EM2_1a(S)
%EM The EM algorithm to estimate parameters
%   f1 = alpha*fc + (1-alpha)*fi1
%   fc is gaussian and fi1 is gumbel
% tollerance = 1e-8;
tollerance = 1e-8;
[N, M] = size(S);

ll = 0;
prev_ll = -1;

s1 = S(1,:);
s1 = s1(s1~=0);

M1 = size(s1,2);

S1_sorted = sort(s1, 'descend');
% S2_sorted = sort(s2, 'descend');
% S_sorted = S;

%prior_thres = 20;
% prior_thres = 28; % pepnovo
prior_thres = 30; % synthetic
nc1 = sum(S1_sorted > prior_thres);
% nc2 = sum(S2_sorted > prior_thres);
% alpha = 0.1;
% alpha = nc1 / M1;
alpha_constraint = 99999; % .10
alpha = min(alpha_constraint,nc1 / M);

[u_c, sigma_c] = normfit(S1_sorted(1:int32(M*alpha)));
gamma_i1 = min(s1);
gamma_i1 = 0;


% start_idx = int32(M * alpha);
% end_idx = M;
% 
% % Ensure indices are within bounds of S1_sorted
% if start_idx < 1
%     start_idx = 1;
% end
% if end_idx > length(S1_sorted)
%     end_idx = length(S1_sorted);
% end
% 
% % Check that start_idx does not exceed end_idx after adjustments
% if start_idx > end_idx
%     error('Adjusted start index exceeds end index. Check the values of M and alpha.');
% end
% 
% % Extract the range and adjust for non-negativity
% selected_data = S1_sorted(start_idx:end_idx) - gamma_i1;
% 
% % Ensure data is non-negative for gamfit
% if any(selected_data < 0)
%     error('Selected data contains negative values after subtracting gamma_i1. Adjust gamma_i1 or ensure S1_sorted values are large enough.');
% end

% Perform the calculation with adjusted, non-negative data
% pa = gamfit(selected_data);


pa = gamfit(S1_sorted(int32(M*alpha):M)-gamma_i1);
a_i1 = pa(1);
b_i1 = 1 / pa(2);

% [a_i1, b_i1] = gumbel_ab(u_i1, sigma_i1);

% figure('Position', [10,10,2000,500]);

stepsize = 0.0001;

% figure;
while abs(ll - prev_ll) > tollerance
    prev_ll = ll;

    ll = func_ll2_1a(s1, alpha, u_c, sigma_c, a_i1, b_i1, gamma_i1);
    disp(ll);
    disp(ll - prev_ll);
    disp([alpha, u_c, sigma_c, a_i1, b_i1, gamma_i1]);

    Rs = rs_gamma(s1, alpha, u_c, sigma_c, a_i1, b_i1, gamma_i1);
    Rsc = Rs;
%     scatter(s1,Rs)
    % plot_dist_gamma_fn(S', alpha, u_c, sigma_c, a_i1, b_i1, gamma_i1);
    % pause(.01);
%     Rsi = 1 - Rs;
    sum_Rsc = sum(Rsc);
%     sum_Rsi = sum(Rsi);
    % alpha_new = sum_Rsc / M;
    % PepNovo: Update alpha
    alpha_new = sum_Rsc / M;
    % alpha_new = min(alpha_new, alpha_constraint); % Enforce prior knowledge constraint (PepNovo sucks)

    u_c_new = (sum( Rsc .* s1 )) / (sum_Rsc);
    sigma_c_new = sqrt(sum( Rsc .* (s1 - u_c_new).^2 ) / (sum_Rsc));
    
%     a_i1_new = b_i1 * (log((M1 - sum_Rsc) / sum( (1 - Rsc) .* exp(-s1 / b_i1) )));
%     b_i1_new = sum( (1 - Rsc) .* (s1 - a_i1_new) .* (1 - exp(-(s1 - a_i1_new)/b_i1)) ) / (M1 - sum_Rsc);
%     
    m1 = sum( (1 - Rsc) .* (s1 - gamma_i1) ) / (M - sum_Rsc);
    m2 = sum( (1 - Rsc) .* (s1 - gamma_i1 - m1).^2 ) / (M - sum_Rsc);
    a_i1_new = m1^2 / m2;
    b_i1_new = m1 / m2;
    
    alpha = alpha_new;
%     gamma_i1 = gamma_i1_new;
    
    u_c = u_c_new;
    sigma_c = sigma_c_new;
    
    a_i1 = a_i1_new;
    b_i1 = b_i1_new;
    
%     [u_i1, sigma_i1] = gumbel_usigma(a_i1, b_i1);
    
end

end
