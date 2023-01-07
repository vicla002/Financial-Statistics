%% 2 
clear all 
clc
T = readtable('BTC-USD.csv');
close = T.Close;
time = T.Date;

%% Box-cox normplot

lambda = bcNormPlot(close)
% Indicates Log-transformation

log_close = log(close);
%% plot data
plot(time, log_close);

%% Differentiate and Check Acf n Pacf
A = [1 -1];
C = [1];
data = filter(A, C, log_close); data = data(length(A):end); %time = time(length(A):end);
plotACFnPACF(data,20,"log closing price for BTC",0.05);

%% Check if data is white
whitenessTest(data,0.05,40);

% Whiteness-test shows data is not white (but very close lol)

%%

ddata = iddata(data);

model_init = idpoly(A, [], C);

% Determine model structure
modelA1 = pem(ddata, model_init);

% Collect residuals
e = resid(ddata, modelA1); e = e(length(A):end);
plotACFnPACF(e.y, 40, "Model residuals", 0.05);
present(modelA1);
var(e.y)


%% Try modeleing
p = 10;
q = 0;
ddata = iddata(data);
model_arma = armax(ddata, [p q]);


% Look into the ACF and the PACF to determine best model structure
e = resid(ddata, model_arma);
plotACFnPACF(e.y, 20, sprintf("ARMA(%d, %d) residuals", p, q), 0.05);

present(model_arma) % Use coefficients as initiation for PEM

var(e.y)

%%
%% PEM ARIMA(10,1,0)
A = [1 1 1 1 1 1 1 1 1 1];
C = 1;

% Differentiation
A_diff = [1 -1];
A_star = conv(A, A_diff);

model_init = idpoly(A_star, [], C);
model_init.Structure.a.Free = [0 ones(1,2) 0 0 0 1 0 0 0 1];

% Determine model structure
modelA2 = pem(ddata, model_init);


% Collect residuals
e = resid(ddata, modelA2); e = e(length(A_star):end);
var(e.y)

plotACFnPACF(e.y, 40, "Model residuals", 0.05);

present(modelA2)

%%

modelA2.A = conv([1 -1], modelA2.A);

%% 1 hour validation data predictions

k = 1;  % Corresponds to one day predictions

% Model A2
% Compute the G and F polynomials.
[F, G] = polydiv(modelA2.C, modelA2.A, k);   

% Form the predicted data for Model A2
log_yhatA2_1d  = filter(G, modelA2.C, log_close);

yhatA2_1d = exp(log_yhatA2_1d);

plot(time, close, "DisplayName", ...
     "Validation data");

hold on
plot(time, yhatA2_1d, "DisplayName", ...
     "Model A2 predictions (1d)");

legend()
axis tight
grid on

%%

% Store the residuals 
ehat_A2 = close - yhatA2_1d;
plotACFnPACF(ehat_A2, 40,"", 0.05);
whitenessTest(ehat_A2,0.05,40);


%% Check whiteness
whitenessTest(modelA2_1h_measured_val_e, 0.05, 50);








%% OK new plan: Lets split up the data to see if we can model better Start with part one

%% Part one
clear all 
clc
T = readtable('BTC-USD.csv');
close = T.Close;
time = T.Date;
plot(close);

%% Split the data in two

close_1 = close(1:length(close)/2);
close_2 = close(length(close)/2:end);
log_close_1 = log(close_1);
log_close_2 = log(close_2);

A = [1 -1];
C = 1;

data_1 = filter(A, C, log_close_1); data_1 = data_1(length(A):end);
data_2 = filter(A, C, log_close_2); data_2 = data_2(length(A):end);

%% Acf and Pacf data 1
plotACFnPACF(data_1,40,"First half",0.05);

%% Acf and Pacf data 2
plotACFnPACF(data_2,40,"Second half", 0.05);

%% Whiteness test
whitenessTest(data_1,0.05,40);
whitenessTest(data_2,0.05,40);

%%



















