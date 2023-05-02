clc; clear;
% Load simulated data
load('data\sim.mat')

%% 
% Box-constraint optimization
ub  =   repmat([ -0.05,  0.09, 40,     5,   30,   30,   -0.01,   8], 1,2);
lb  =   repmat([ -0.25, 0.005,  1,   0.05, 0.005,   0,   -0.1,  -5], 1,2);

% Sparse-gridd quadrature rule for numerical integration
[ Int, WTS, PTS, INTCLS ] = fwtpts( 2, 11, 'Norm');
[ Int, WTS4, PTS4, INTCLS ] = fwtpts( 4, 11, 'Norm');

crit = @(theta)(mSVhatHJ_crit_inst4(theta,...
    PTS, WTS', mY, mOptPriceIV1, mOptPriceIV2, mV(2:end,:), mK1, mK2, vTau, dt, r, 50));

% Non-linear constraint
constr= @(theta)(mSVhatHJ_fmin_constr(theta));
optSearch = optimset('Display','iter', 'PlotFcns', @optimplotfval)

% Optimization
[vTheta1, fval] = fminsearchcon(crit,vStart,lb,ub,[],[],constr,optSearch);

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

% Figure
mIntens = mSVhatHJ_ImpIntens(mOptPriceIV1, mOptPriceIV2, mY, mV, mK1, mK2, vTau, r, mTheta);

subplot(2,1,1); plot(mIntens(:,1));
title("Implied intensity - stock 1")
subplot(2,1,2); plot(mIntens(:,2));
title("Implied intensity - stock 2")


%%
save('data\sim.mat')


