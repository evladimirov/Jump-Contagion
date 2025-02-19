function [mSolODE] = mSVhatHJ_ODE_cos(T, r, mParam)
%   ODE solver for 2SVHJ model for all arguments from COS method option pricing
% 
%       Inputs:
%           T       double, time-to-maturity, in YEARs
%           r       double, risk-free rate
%           vParam  vector of the following parameters:
%                     sigma_v, kappa_v, vbar, rho, muj_q, sigmaj, kappa_l, lbar, delta
%       Return:
%           mSol    Nx4 matrix, solutions of ODE for SVHJ model for all
%                   vKu arguments from Carr-Madan
%
%   author: Evgenii Vladimirov
%   date:   19.03.2019
%   last version: 30.04.2019
%
%%
    
    [muj_q1, sigmaj1, kl1, lbar1, delta1, mut_delta1] = deal(mParam(1,1),mParam(1,2),mParam(1,3), mParam(1,4), mParam(1,5), mParam(1,6));
    [    ~ ,      ~ , kl2, lbar2, delta2, mut_delta2] = deal(mParam(2,1),mParam(2,2),mParam(2,3), mParam(2,4), mParam(2,5), mParam(2,6));
    
    % Define CCF system matrices for the model
    K0 = [ 0;
         kl1*lbar1;
         kl2*lbar2];
    K1 = [ 0,   -(exp(muj_q1 + 0.5*sigmaj1^2)-1),   0;
           0,      -kl1,                            0;
           0,       0,                             -kl2];

    H0 = zeros(3,3);
    H1 = zeros(3,3,3);

    L0 = [0, 0];
    L1 = [0, 0;
          1, 0;
          0, 1];
    
    JT = @(vBeta)([exp(muj_q1*vBeta(1) +0.5*sigmaj1^2*vBeta(1)^2 +vBeta(2)*delta1 + vBeta(3)*mut_delta2)-1; 
                   exp(vBeta(2)*mut_delta1 + vBeta(3)*delta2)-1]);

    % COS based routine
    a  = -5.0;               % Integration lower bound
    b  =  5.0;               % Integration upper bound
    N  = 1024;               % Integration bound for COS pricing part
    k = cumsum(ones(N,1))-1; % Integration grid values  

    mSolODE = zeros(N,4,length(T));
    parfor i=1:N
        mSolODE(i,:,:) = affineODE([k(i)*pi/(b-a);0;0], T, K0, K1, H0, H1, L0, L1, JT);
    end     
end

