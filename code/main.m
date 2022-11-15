clc; clear;
% Load simulated data
load('share\data\sim.mat')

%% 
% Box-constraint optimization
ub  =   repmat([ -0.05,  0.09, 40,     5,   30,   30,   -0.01,   8], 1,2);
lb  =   repmat([ -0.25, 0.005,  1,   0.05, 0.005,   0,   -0.1,  -5], 1,2);

% Sparse-gridd quadrature rule for numerical integration
[ Int, WTS, PTS, INTCLS ] = fwtpts( 2, 11, 'Norm');

crit = @(theta)(mSVhatHJ_crit_inst4(theta,...
    PTS, WTS', mY, mOptPriceIV1(:,idK,idTau), mOptPriceIV2(:,idK,idTau), mV(2:end,:), ...
    mK1(:,idK), mK2(:,idK), vTau(idTau),dt,r, 50));

% Non-linear constraint
constr= @(theta)(mSVhatHJ_fmin_constr(theta));
optSearch = optimset('Display','iter', 'PlotFcns', @optimplotfval)

% Optimization
[vTheta, fval] = fminsearchcon(crit,vStart,lb,ub,[],[],constr,optSearch)

display(vTheta)

