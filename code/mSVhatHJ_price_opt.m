function [mOptPrice1, mOptPrice2] = mSVhatHJ_price_opt(mY, mV, mLambda, r, mParam, T, M, bBSIV)
%  Price option panel from 2SVhatHJ model
%
%   Used to produce Table R.2 in the replication pdf (analogous to Table 4 in the main text)
%
%       Inputs:
%           iN          integer, number of draws
%           mY          iNx2 matrix of log stock prices
%           mV          iNx2 matrix of volatilities
%           mLambda     iNx2 matrix of intensities
%           r           double, risk-free rate
%           mParam      matrix of parameters
%           T           column-vector, time-to-maturity
%           M           vector, moneyness: stocks by rows
%           bBSIV       boolean, if 1 -> option prices in currency
%                                if 0 -> option prices in BSIV terms
%
%       Return:
%           mOptPrice1  (iN+1)xK matrix of option prices
%           mOptPrice2  (iN+1)xK matrix of option prices
%
%   author: Evgenii Vladimirov
%   date:   19.03.2019
%
%%

    mSolODE1 = mSVhatHJ_ODE_cos(T, r, mParam);
    mSolODE2 = mSVhatHJ_ODE_cos(T, r, flip(mParam));
    
    iK = length(M);
    iN = size(mY,1);
    mOptPrice1 = zeros(iN,iK);
    mOptPrice2 = mOptPrice1;

    u_tmp = mSolODE1(:,2);

    SolODE1_withVol = repmat(mSolODE1, 1,1,iN);
    SolODE2_withVol = repmat(mSolODE2, 1,1,iN);

    SolODE1_withVol(:,1,:) = squeeze(SolODE1_withVol(:,1,:)) +  0.5*mV(:,1)'.*( -u_tmp + u_tmp.^2 )*T;
    SolODE2_withVol(:,1,:) = squeeze(SolODE2_withVol(:,1,:)) +  0.5*mV(:,2)'.*( -u_tmp + u_tmp.^2 )*T;

    parfor kk = 1:iK
        for i = 1:iN
            mOptPrice1(i,kk) = Opt_price_cos([mY(i,1);mLambda(i,1);mLambda(i,2)], M(kk)*exp(mY(i,1)), T, r, [], 1, SolODE1_withVol(:,:,i), bBSIV);
            mOptPrice2(i,kk) = Opt_price_cos([mY(i,2);mLambda(i,2);mLambda(i,1)], M(kk)*exp(mY(i,2)), T, r, [], 1, SolODE2_withVol(:,:,i), bBSIV);
        end
    end
end

