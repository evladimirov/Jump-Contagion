function [mSolODE] = MSVHJ_ODE_cos(T, r, mParam)
%  ODE solver for 2SVHJ model for all arguments from COS option pricing
% 
%  Used to produce Table R.5 in the replication pdf (Table 2 and Table C.1 in the main text)
%
%       Inputs:
%           T       double, time-to-maturity, in YEARs
%           r       double, risk-free rate
%           mParam  matrix of the parameters:
%
%       Return:
%           mSol    Nx4 matrix, solutions of ODE for SVHJ model for all
%                   arguments from the COS method
%
%   author: Evgenii Vladimirov
%   date:   19.03.2019
%   last version: 30.04.2019
%
%%
    [sigma_v, kv, vbar, rho] = deal(mParam(1,1),mParam(1,2), mParam(1,3), mParam(1,4));
    [muj_q, sigmaj, kl1, lbar1, delta1, mut_delta1] = deal(mParam(1,5),mParam(1,6),mParam(1,7), mParam(1,8), mParam(1,9), mParam(1,10));
    [kl2, lbar2, delta2, mut_delta2] = deal(mParam(2,7), mParam(2,8), mParam(2,9), mParam(2,10));
    
    %Define CCF system matrices
    K0 = [ r;
         kv*vbar;
         kl1*lbar1;
         kl2*lbar2];
    K1 = [ 0,   -0.5,   -(exp(muj_q + 0.5*sigmaj^2)-1),   0;
           0,   -kv,      0,                              0;
           0,    0,      -kl1,                            0;
           0,    0,       0,                             -kl2];

    H0 = zeros(4,4);
    H1 = zeros(4,4,4);
    H1(:,:,2) = [ 1,          sigma_v*rho,  0,  0; 
                 sigma_v*rho, sigma_v^2,    0,  0;
                 0,           0,            0,  0;
                 0,           0,            0,  0];
    L0 = [0, 0];
    L1 = [0, 0;
          0, 0;
          1, 0;
          0, 1];
    
    JT = @(vBeta)([exp(muj_q*vBeta(1) +0.5*sigmaj^2*vBeta(1)^2 +vBeta(3)*delta1 + vBeta(4)*mut_delta2)-1; 
                   exp(vBeta(3)*mut_delta1 + vBeta(4)*delta2)-1]);

    %COS based routine
    a  = -5.0;   %integration lower bound
    b  =  5.0;   %integration upper bound
    N	= 1024;  %Integration bound for COS pricing part
    k = cumsum(ones(N,1))-1; %Integration grid values  

    mSolODE = zeros(N,5, length(T));
    parfor i=1:N
        mSolODE(i,:,:) = affineODE([k(i)*pi/(b-a);0;0;0], T, K0, K1, H0, H1, L0, L1, JT);
      %the last arguments passed to ODE set = 0 as we need CF of log(S) only, not a joint CF
    end
               
end

