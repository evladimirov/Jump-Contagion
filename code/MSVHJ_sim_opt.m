function [mOptPrice1, mOptPrice2] = MSVHJ_sim_opt(iN, mY, mV, mLambda, r, mParam, T, M, bBSIV)
%  Simulate option panel from 2SVHJ model
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
%           mOptPrice1  (iN+1)xKxT 3d matrix of option prices
%           mOptPrice2  (iN+1)xKxT 3d matrix of option prices
%
%   author: Evgenii Vladimirov
%   date:   19.03.2019
%
%%

    mSolODE1 = MSVHJ_ODE_cos(T, r, mParam);
    mSolODE2 = MSVHJ_ODE_cos(T, r, flip(mParam));
    
    iK = length(M);
    iT = length(T);
    mOptPrice1 = zeros(iN+1,iK,iT);
    mOptPrice2 = mOptPrice1;
    parfor kk = 1:iK
        for tt = 1:iT
            for i=1:(iN+1)
                mOptPrice1(i,kk,tt) = Opt_price_cos([mY(i,1);mV(i,1);mLambda(i,1);mLambda(i,2)], M(kk)*exp(mY(i,1)), T(tt), r, [], 1, mSolODE1(:,:,tt), bBSIV);
                mOptPrice2(i,kk,tt) = Opt_price_cos([mY(i,2);mV(i,2);mLambda(i,2);mLambda(i,1)], M(kk)*exp(mY(i,2)), T(tt), r, [], 1, mSolODE2(:,:,tt), bBSIV);
            end
        end
    end
end

