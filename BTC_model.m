clear all
clear var

%% 2 

%adjClose inkluderar stock splits, dividends etc medan Closed bara är pris
T = readtable('BTC-USD.csv');
average = zeros(1949,1);
high = T.High;
low = T.Low;
open = T.Open;
close = T.Close;
for i=1:1949
    average(i,1) = (high(i)+low(i))/2;
end

for i=1:1949
    averageDay(i,1) = (open(i)+close(i))/2;
end

%plot(averageDay-average)
returnsClose = zeros(1,1948);
for i=2:1949
    returnsClose(i)=close(i)-close(i-1);
end

logReturnsClose = zeros(1948,1);
for i=1:1948
    logReturnsClose(i)=log(returnsClose(i));
end
plot(returnsClose)

omega = 1;
alpha = 1;
beta = 1;
[params_garch, lOut_garch, CovM_garch] = MLmax(@lnL_garch, [omega alpha beta mean(returnsClose)], returnsClose);
conf_int_omega = [params_garch(1) - 1.96 * sqrt(CovM_garch(1, 1)) params_garch(1) + 1.96 * sqrt(CovM_garch(1, 1))]
conf_int_alpha = [params_garch(2) - 1.96 * sqrt(CovM_garch(2, 2)) params_garch(2) + 1.96 * sqrt(CovM_garch(2, 2))]
conf_int_beta = [params_garch(3) - 1.96 * sqrt(CovM_garch(3, 3)) params_garch(3) + 1.96 * sqrt(CovM_garch(3, 3))]
conf_int_my = [params_garch(4) - 1.96 * sqrt(CovM_garch(4, 4)) params_garch(4) + 1.96 * sqrt(CovM_garch(4, 4))]
 
omega_1 = params_garch(1);
alpha_1 = params_garch(2);
beta_1 = params_garch(3);
mu_1 = params_garch(4);

figure
plot(returnsClose)
plot(logReturnsClose)