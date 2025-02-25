function [mCF] = SVhatHJ_ccf2(s, mX, vVspot, dt, vParam, c)
%     CCF for SVhatHJ model
% 
%     Inputs:
%         s         2xP matrix of stacked arguments
%         mX        Nx2 matrix of states, i.e. stacked state vectors
%                                       (log(S_t), lambda_t)
%         vVspot    Nx1 vector of spot volatility estimates
%         dt        time discretisation
%         vParam    vector of the following parameters:
%                    muj_q, sigmaj, kappa_l, lbar, delta, muj, eta
%
%     Output:
%         mCF       matrix of CCF values for SVhatHJ model
%
%   author: Evgenii Vladimirov
%   date:   30.03.2019 
%
%%
    [muj_q, sigmaj, kappa_l, lbar, delta, muj, eta] = deal(vParam(1),vParam(2),vParam(3), vParam(4), vParam(5), vParam(6), vParam(7));

    %Define CCF system matrices
    K0 = [ 0;
         kappa_l*lbar];
    K1 = [ 0,   -c*(exp(muj_q + 0.5*sigmaj^2)-1);
           0,   -kappa_l];

    H0 = zeros(2,2);
    H1 = zeros(2,2,2);
    L0 = 0;
    L1 = [0; 1];
    JT = @(vBeta)(exp(muj*vBeta(1)*c + 0.5*(sigmaj)^2*vBeta(1)^2*c^2 + vBeta(2)*delta)-1); 

    iP = size(s,2);
    mSolODE = zeros(3,iP);
    
    parfor j =1:iP
        mSolODE(:,j) = transpose(affineODE(s(:,j), dt, K0, K1, H0, H1, L0, L1, JT));
    end
    
    mCF = exp((vVspot(1:end-1)*((eta-0.5)*1i*s(1,:)*c - 0.5*c^2*s(1,:).^2))*dt + mSolODE(1,:) + mX(1:end-1,2)*mSolODE(3,:));

end

