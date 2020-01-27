function [u,sigma,lambda] = sn_para_est(X)
%sn_para_est estimate the parameters of the skew normal with method of
%moments
%   Detailed explanation goes here

m1 = mean(X);
m2 = var(X);
m3 = skewness(X);
a1 = sqrt(2/pi);
b1 = (4/pi - 1) / a1;

delta = sign(m3) ./ sqrt(a1^2 + m2 .* (b1 / abs(m3)).^(2/3));

sigma = sqrt(m2 ./ (1 - a1^2 * delta.^2));

u = m1 - a1*delta*sigma;

lambda = m3;

end

