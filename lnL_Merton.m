function prob = lnL_Merton(thetas, indata)

% Number of jumps to check
n_max = 360;
n_poisson_part= zeros(n_max, 1);
sigma_2_vec = zeros(n_max, 1);
mu_vec = zeros(n_max, 1);

% Lambda, sigma, and delta have to be larger than 0.
gamma = thetas(1);
sigma = abs(thetas(2));
lambda = abs(thetas(3)); 
mu = thetas(4);
delta = abs(thetas(5));
delta_t = 1/12;

for n = 0:n_max-1 
    sigma_2_vec(n + 1, 1) = sigma^2*delta_t + n*delta^2;
    mu_vec(n + 1, 1) = gamma*delta_t + n*mu;
    n_poisson_part(n + 1, 1) = ((lambda*delta_t)^n/factorial(n)) * exp(-lambda*delta_t);
end

prob = normpdf(indata, mu_vec', sqrt(sigma_2_vec)')*n_poisson_part;
%prob = 1./sqrt(2.*pi.*sigma_2_vec).*exp((indata-mu)./sigma_2_vec).*n_poisson_part;
prob = log(prob);



