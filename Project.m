%% Task 1a
clc 
clear all
load MertonDataA.mat
gamma = 1;
sigma = 1;
lambda = 1;
my = 1;
delta = 1;
theta = [gamma, sigma, lambda, my, delta];

[paramsOut_Merton, lOut_Merton, CovM_Merton] = MLmax(@lnL_Merton, theta, XA);

gamma_1 = paramsOut_Merton(1);
sigma_1 = paramsOut_Merton(2);
lambda_1 = paramsOut_Merton(3);
my_1 = paramsOut_Merton(4);
delta_1 = paramsOut_Merton(5);


%% GMM
medel = mean(XA);
XA = XA-mean(XA);  % vad ska egentligen yN vara

beta0 = 1;
beta1 = 1;
a = 1;
b = 1;
c = 1;
theta0 = [beta0, beta1, a, b, c];
iterations = 1;
W = eye(4,4);
GMM1 = @(theta)GMM2(theta0,XA,W)
theta_hat = fminunc(GMM1,theta0);

for i = 1:iterations
    [m, w_new] = GMM2(theta0,XA,W)
    w = w_new;
    GMM1 = @(theta)GMM2(theta,XA,W)
    theta_hat = fminunc(GMM1,theta0)
end



%% Koden som ish funkar
clear all 
clc
close all
T = readtable('BTC-USD.csv');
[h,pValue,stat,cValue] = adftest(T.Close)
%%
figure
plot(T.Close)
saveas(gcf,'figs/priceProcess.png')
returns = zeros(1,1948);

for i = 2:1949
    returns(i) = (T.Close(i)-T.Close(i-1))/T.Close(i-1);
end
figure
plot(returns)
saveas(gcf, 'figs/dailyReturns.png')
omega = 0.1;
alpha = 0.25;
beta = 0.6; 
my = mean(returns);

[params_garch, lOut_garch, CovM_garch] = MLmax(@lnL_garch, [omega alpha beta my], returns);
conf_int_omega = [params_garch(1) - 1.96 * sqrt(CovM_garch(1, 1)) params_garch(1) + 1.96 * sqrt(CovM_garch(1, 1))]
conf_int_alpha = [params_garch(2) - 1.96 * sqrt(CovM_garch(2, 2)) params_garch(2) + 1.96 * sqrt(CovM_garch(2, 2))]
conf_int_beta = [params_garch(3) - 1.96 * sqrt(CovM_garch(3, 3)) params_garch(3) + 1.96 * sqrt(CovM_garch(3, 3))]
conf_int_my = [params_garch(4) - 1.96 * sqrt(CovM_garch(4, 4)) params_garch(4) + 1.96 * sqrt(CovM_garch(4, 4))]
 
omega_1 = params_garch(1);
alpha_1 = params_garch(2);
beta_1 = params_garch(3);
mu_1 = params_garch(4);

sigmas = zeros(length(returns),1);
for i =2:length(sigmas)
    sigmas(i) = omega_1 + alpha_1 * (returns(i - 1)-mu_1)^2 + beta_1 * sigmas(i - 1);
end
plot(sigmas.^2)
figure
plot(returns.^2)




%% ARMA
clear all 
clc
T = readtable('BTC-USD.csv');
autocorr(T.Close)
returns = zeros(1,1948);
plot(T.Close)
for i = 2:1949
    returns(i) = (T.Close(i)-T.Close(i-1))/T.Close(i-1);
end
figure
autocorr(returns)
figure
plot(T.Close)

detrender = [1 -1];
T.close_detrend = filter(detrender, 1,T.Close);
hold on
plot(T.close_detrend)
legend(["Original", "Detrended"])

%%
autocorr(T.close_detrend)
figure
autocorr(returns)

p = ones(1,9);
q = ones(1,3);
[params_ARMA, lOut_ARMA, CovM_ARMA] = MLmax(@lnL_ARMA, [p, q], T.close_detrend');
%%














%% MATLABS GARCH SIMULERING
clear all 
clc
T = readtable('BTC-USD.csv');
close = T.Close;
time = T.Date;
tt = timetable(time,close);


%%

model = GJR_pricesLogReturn
V0 = infer(model,pricesLogReturn);
[V,Y] = simulate(model,365,'NumPaths',10000,'E0', pricesLogReturn, 'V0', V0);

%%
figure
subplot(2,1,1)
plot(V(:,1:10))
title('Simulated Conditional Variances')
subplot(2,1,2)
plot(Y(:,1:10))
title('Simulated Returns')
%%
condVol = sqrt(V0);
plot(pricesLogReturn); hold on;
plot(condVol); hold off;

%% Koden som ish funkar
clear all 
clc
T = readtable('BTC-USD.csv');
returns = zeros(1,1948);
for i = 2:1949
    returns(i) = (T.Close(i)-T.Close(i-1))/T.Close(i-1);
end
%plot(returns)
omega = 0.1;
alpha = 0.25;
beta = 0.6; 
my = mean(returns);

[params_garch, lOut_garch, CovM_garch] = MLmax(@Copy_of_lnL_garch, [omega alpha beta my], returns);
conf_int_omega = [params_garch(1) - 1.96 * sqrt(CovM_garch(1, 1)) params_garch(1) + 1.96 * sqrt(CovM_garch(1, 1))]
conf_int_alpha = [params_garch(2) - 1.96 * sqrt(CovM_garch(2, 2)) params_garch(2) + 1.96 * sqrt(CovM_garch(2, 2))]
conf_int_beta = [params_garch(3) - 1.96 * sqrt(CovM_garch(3, 3)) params_garch(3) + 1.96 * sqrt(CovM_garch(3, 3))]
conf_int_my = [params_garch(4) - 1.96 * sqrt(CovM_garch(4, 4)) params_garch(4) + 1.96 * sqrt(CovM_garch(4, 4))]
 
omega_1 = params_garch(1);
alpha_1 = params_garch(2);
beta_1 = params_garch(3);
mu_1 = params_garch(4);

sigmas = zeros(length(returns),1);
for i =2:length(sigmas)
    sigmas(i) = omega_1 + alpha_1 * (returns(i - 1)-mu_1)^2 + beta_1 * sigmas(i - 1);
end

figure
plot(sigmas)
hold on
plot(returns.^2)
saveas(gcf,'figs/test.png')

%% Testar Axel Sjöberg version 
clear all 
clc
close all
T = readtable('BTC-USD.csv');
returns = zeros(1,1948);
for i = 2:1949
    returns(i) = log((T.Close(i)/(T.Close(i-1))));
end
returns = returns-mean(returns);

omega = 0.1;
alpha = 0.25;
beta = 0.6; 
my = mean(returns);

[params_garch, lOut_garch, CovM_garch] = MLmax(@lnL_garch, [omega alpha beta my], returns);
params_garch = abs(params_garch)
conf_int_omega = [params_garch(1) - 1.96 * sqrt(CovM_garch(1, 1)) params_garch(1) + 1.96 * sqrt(CovM_garch(1, 1))]
conf_int_alpha = [params_garch(2) - 1.96 * sqrt(CovM_garch(2, 2)) params_garch(2) + 1.96 * sqrt(CovM_garch(2, 2))]
conf_int_beta = [params_garch(3) - 1.96 * sqrt(CovM_garch(3, 3)) params_garch(3) + 1.96 * sqrt(CovM_garch(3, 3))]
conf_int_my = [params_garch(4) - 1.96 * sqrt(CovM_garch(4, 4)) params_garch(4) + 1.96 * sqrt(CovM_garch(4, 4))]
 
omega_1 = abs(params_garch(1));
alpha_1 = abs(params_garch(2));
beta_1 = abs(params_garch(3));
mu_1 = abs(params_garch(4));

sigmas = zeros(length(returns),1);
sigmas(1) = var(returns);
for i =2:length(sigmas)
    sigmas(i) = omega_1 + alpha_1 * (returns(i - 1)-mu_1)^2 + beta_1 * sigmas(i - 1);
end
figure
plot(sigmas)
hold on
plot(returns.^2)

%fan det funkar ju TYP med abs! 


%% Testar garch 2,2 
clear all 
clc
close all
T = readtable('BTC-USD.csv');
returns = zeros(1,1948);
for i = 2:1949
    returns(i) = log((T.Close(i)/(T.Close(i-1))));
end
returns = returns-mean(returns);

omega = 0.1;
alpha = 0.25;
beta = 0.6; 
alpha2 = 0.25;
beta2 = 0.6;
my = mean(returns);

[params_garch, lOut_garch, CovM_garch] = MLmax(@lnL_garch12, [omega alpha beta my alpha2 beta2], returns);
conf_int_omega = [params_garch(1) - 1.96 * sqrt(CovM_garch(1, 1)) params_garch(1) + 1.96 * sqrt(CovM_garch(1, 1))]
conf_int_alpha = [params_garch(2) - 1.96 * sqrt(CovM_garch(2, 2)) params_garch(2) + 1.96 * sqrt(CovM_garch(2, 2))]
conf_int_beta = [params_garch(3) - 1.96 * sqrt(CovM_garch(3, 3)) params_garch(3) + 1.96 * sqrt(CovM_garch(3, 3))]
conf_int_my = [params_garch(4) - 1.96 * sqrt(CovM_garch(4, 4)) params_garch(4) + 1.96 * sqrt(CovM_garch(4, 4))]
 
omega_1 = abs(params_garch(1));
alpha_1 = abs(params_garch(2));
beta_1 = abs(params_garch(3));
mu_1 = abs(params_garch(4));

sigmas = zeros(length(returns),1);
for i =2:length(sigmas)
    sigmas(i) = omega_1 + alpha_1 * (returns(i - 1)-mu_1)^2 + beta_1 * sigmas(i - 1);
end
figure
plot(sigmas)

