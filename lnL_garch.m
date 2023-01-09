function l=lnL_garch(params,r)
% Output must be a column vector

omega = abs(params(1));
alpha = abs(params(2));
beta = abs(params(3));
mu = abs(params(4));
N = length(r);

%create sigma vector 
sigma2 = zeros(1,N);
sigma2(1) = var(r);

for i = 2:N
    sigma2(i) = omega + alpha * (r(i - 1)-mu)^2 + beta * sigma2(i - 1);
end

l=-1/2*log(2*pi*sigma2)-1/2*((r-mu).^2)./sigma2;
