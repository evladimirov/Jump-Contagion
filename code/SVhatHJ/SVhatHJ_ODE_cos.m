function [mSolODE] = SVhatHJ_ODE_cos(T, r, vParam)
%  ODE solver for modified SVhatHJ model for all arguments from Fang-Oosterlee COS option pricing
%  CCF matrices are modified s.t. alpha has to be adjusted for full solution
%
%  Used to produce Table R.3 in the replication pdf (analogous to Table 5 in the main text)
%
%       Inputs:
%           T       column-vector, time-to-maturity, in YEARs
%           r       double, risk-free rate
%           vParam  vector of the following parameters:
%                      muj_q, sigmaj, kappa_l, lbar, delta
%       Return:
%           mSol    Nx4 matrix, solutions of ODE for SVhatHJ model for all
%                   arguments from the COS method
%
%   author: Evgenii Vladimirov
%   date:   29.03.2019
%   last version: 27.04.2019
%
%%
    [muj_q, sigmaj, kappa_l, lbar, delta] = deal(vParam(1),vParam(2),vParam(3), vParam(4), vParam(5));
    
    %Define CCF system matrices for the model
    K0 = [ 0;
         kappa_l*lbar];
    K1 = [ 0,   -(exp(muj_q + 0.5*sigmaj^2)-1);
           0,   -kappa_l];

    H0 = zeros(2);
    H1 = zeros(2,2,2);

    L0 = 0;
    L1 = [0;1];
    JT = @(vBeta)(exp(muj_q*vBeta(1) +0.5*sigmaj^2*vBeta(1)^2 +vBeta(2)*delta)-1); 

    %Cosine based routine
    a  = -5.0;                   % integration lower bound
    b  =  5.0;                   % integration upper bound
    N_cos 	= 1024;              % integration bound for COS pricing part
    k = cumsum(ones(N_cos,1))-1; % integration grid values  

    mSolODE = zeros(N_cos,3,length(T));
    parfor i=1:N_cos
        mSolODE(i,:,:) = affineODE([k(i)*pi/(b-a);0], T, K0, K1, H0, H1, L0, L1, JT);
        %2nd argument passed to ODE set = 0 as we need CF of log(S) only
    end
               
end

