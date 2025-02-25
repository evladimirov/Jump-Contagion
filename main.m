 clc; clear;
% Load simulated data
load('data/sim.mat')

%% Estimation. Bivariate model

% Box-constraint optimization
ub  =   repmat([ -0.05,  0.09, 40,     5,   30,   30,   -0.01,   8], 1,2);
lb  =   repmat([ -0.25, 0.005,  1,   0.05, 0.005,   0,   -0.1,  -5], 1,2);

% Sparse-grid quadrature rule for numerical integration
[ ~, WTS, PTS, ~ ] = fwtpts( 2, 11, 'Norm');
[ Int, WTS4, PTS4, INTCLS ] = fwtpts( 4, 11, 'Norm');

crit = @(theta)(mSVhatHJ_crit_inst4(theta,...
    PTS, WTS', mY, mOptPriceIV1, mOptPriceIV2, mV(2:end,:), mK1, mK2, vTau, dt, r, 50));

% Non-linear constraint
constr= @(theta)(mSVhatHJ_fmin_constr(theta));
optSearch = optimset('Display','iter', 'PlotFcns', @optimplotfval);

% Optimization
[vTheta, fval] = fminsearchcon(crit,vStart,lb,ub,[],[],constr,optSearch);

mTheta = [vTheta(1:end/2); vTheta(end/2+1:end)];
disp('Estimated parameters:')
disp(mTheta)

% Standard Errors
[mA, mB] = mSVhatHJ_std4(vTheta, PTS4, WTS4, mY, mOptPriceIV1, mOptPriceIV2,...
        mV, mK1, mK2, vTau, dt, r, 50);
vStd = real(sqrt(diag(pinv(mA)*mB*pinv(mA)/iN)))';

mStd = [vStd(1:end/2); vStd(end/2+1:end)];
disp('Standard errors:')
disp(mStd)

%% Figure of jump intensities

mIntens = mSVhatHJ_ImpIntens(mOptPriceIV1, mOptPriceIV2, mY, mV, mK1, mK2, vTau, r, mTheta);

subplot(2,1,1); plot(mIntens(:,1));
title("Implied intensity - stock 1")
subplot(2,1,2); plot(mIntens(:,2));
title("Implied intensity - stock 2")


%% Option prices: fit given simulated data

mParam = [vParam(1:end/2); vParam(end/2+1:end)];
vM = [0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1];

[mOptPriceIV1_hat, mOptPriceIV2_hat] = mSVhatHJ_price_opt(mY, mV, mIntens, r, mTheta, vTau, vM, 1);

rmse1 = [100*sqrt(mean((mOptPriceIV1_hat - mOptPriceIV1).^2,1)), 100*sqrt(mean((mOptPriceIV1_hat - mOptPriceIV1).^2, 'all'))];
rmse2 = [100*sqrt(mean((mOptPriceIV2_hat - mOptPriceIV2).^2,1)), 100*sqrt(mean((mOptPriceIV2_hat - mOptPriceIV2).^2, 'all'))];

VarNames= {'k=0.8','k=0.85','k=0.9','k=0.95','k=1.0','k=1.05','k=1.1','total'};
VarNames = matlab.lang.makeValidName(VarNames);
tab_options = array2table(round([rmse1; rmse2],2),'RowNames',{'stock-1','stock-2'},...
        'VariableNames',VarNames);

disp('---------------------------------------------------------------------------');
disp('          Table: Option prices: fit using simulated data)')
disp(' ')
disp(tab_options)

%% Univariate results given simulated data

ub1 = [ -0.05,  0.09, 40,     5,     30,    -0.01,   8];
lb1 = [ -0.25, 0.005,  1,   0.05, 0.005,    -0.1,  1.5];

A1  =   [ 0,    0,     -1,      0,     1,     0   , 0  ];
b1   =   -eps;
optSearch = optimset('Display','iter', 'PlotFcns',@optimplotfval);

% no scaling with backing-out only intensities
crit_univ1 = @(theta)(SVhatHJ_crit_inst(theta, PTS, WTS', mY(:,1), mOptPriceIV1, mV(:,1), mK1, vTau, dt, r, 50));

[vTheta1_univ,fval1_univ] = fminsearchcon(crit_univ1, vStart([1:5,7:8]),lb1,ub1,A1,b1,[],optSearch);

disp('Estimated parameters for univariate model:')
disp(vTheta1_univ)

%% Replicate results from the paper given the empirical estimates

% Table 6
res_tab6 = Table6();

% Figure 3
Figure3();

%% Simulation results: Table 2 and Table C.1
% CAUTION! Simulation takes long time (a few days) to run!

MCTableRes = mSVhatHJ_monte_carlo();


%%
save('data/sim.mat')


