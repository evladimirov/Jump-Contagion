function [mY_hf, mY, mLambda, mCount] = mSVhatHJ_sim(iN, dt, iM, mX0, r, mParam, corr, v1, v2)
%   Simulate state data from the 2SVhatHJ model with constant volatilities 
%
%       Inputs:
%           iN      integer, number of observations
%           dt      double, time discretisation (dt = iT/iN) between two
%                           ovservations
%           iM      integer, number of draws in between two data points
%           mX0     2xK matrix, initial state vector
%           r       double, risk-free rate
%           mParam  KxP matrix of parameters for each of the K stocks
%           corr    double, instanteneous correlation between Brownians in
%                           2 stocks
%           v1      double, variance of stock 1
%           v2      double, variance of stock 2
%
%       Return:
%           mY_hf       (iN+1)*100x2 matrix of HF stock prices
%           mY          (iN+1)x2 matrix of simulated stock path
%           mLambda     (iN+1)x2 matrix of intensities
%
%   author: Evgenii Vladimirov
%   date:   12.10.2020
%
%% 
    % increase the simulation frequency for higher precision
    iP = iM*iN;
    dtt = dt/iM;
    
    mdW = sqrt(dtt)*randn(iP,2);                % components of Brownian motions
    vB = corr*mdW(:,1)+sqrt(1-corr^2)*mdW(:,2); % diffusive component in the second process
    
    mU = rand(iP,2);
     
    mY = [mX0(1,:); zeros(iP,2)];
    mLambda = [mX0(2,:); zeros(iP,2)];
    mCount = [0, 0];

    [muj_q1, sigmaj1, kl1, lbar1, delta1, mut_delta1, muj1, eta1] = deal(mParam(1,1),mParam(1,2),mParam(1,3), mParam(1,4),...
        mParam(1,5), mParam(1,6), mParam(1,7), mParam(1,8));

    [muj_q2, sigmaj2, kl2, lbar2, delta2, mut_delta2, muj2, eta2] = deal(mParam(2,1),mParam(2,2),mParam(2,3), mParam(2,4),...
        mParam(2,5), mParam(2,6), mParam(2,7), mParam(2,8));
    
    % Simulate the log price and jump intensity using the Eurler discretization scheme
    for i=2:(iP+1)

        mLambda(i,1) = mLambda(i-1,1) + kl1*(lbar1 - mLambda(i-1,1))*dtt + delta1*(mLambda(i-1,1)*dtt > mU(i-1,1)) + mut_delta1*(mLambda(i-1,2)*dtt > mU(i-1,2));
        mLambda(i,2) = mLambda(i-1,2) + kl2*(lbar2 - mLambda(i-1,2))*dtt + mut_delta2*(mLambda(i-1,1)*dtt > mU(i-1,1)) + delta2*(mLambda(i-1,2)*dtt > mU(i-1,2));

        mY(i,1) = mY(i-1,1) + ((eta1-0.5)*v1-(exp(muj_q1+0.5*sigmaj1^2)-1)*mLambda(i-1,1))*dtt ...
            + sqrt(v1)*mdW(i-1,1) + (mLambda(i-1,1)*dtt > mU(i-1,1))*(muj1+sigmaj1*randn(1));
        mY(i,2) = mY(i-1,2) + ((eta2-0.5)*v2-(exp(muj_q2+0.5*sigmaj2^2)-1)*mLambda(i-1,2))*dtt ...
            + sqrt(v2)*vB(i-1) + (mLambda(i-1,2)*dtt > mU(i-1,2))*(muj2+sigmaj2*randn(1));

        mCount(i,1) = mCount(i-1,1) + (mLambda(i-1,1)*dtt > mU(i-1,1));
        mCount(i,2) = mCount(i-1,2) + (mLambda(i-1,2)*dtt > mU(i-1,2));
    end

    % Subsample for the desired frequency
    mY_hf = mY;
    mY = mY(1:iM:iP+1,:);
    mLambda= mLambda(1:iM:iP+1,:);
    mCount = mCount(1:iM:iP+1,:);
end

