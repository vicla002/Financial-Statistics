%% 2 
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