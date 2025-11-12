function [out] = SVhatHJ_int_inst(s, WTS, mX, vVspot, dt, vParam,  c)
%     Integrand of the criterion function for the 1st step C-GMM for SVhatHJ model
% 
%  Used to produce Table R.3 in the replication pdf (analogous to Table 5 in the main text)
%
%     Inputs:
%         s         2xP matrix of stacked arguments
%         mX        Nx2 matrix of states, i.e. stacked state vectors
%                                       (log(S_t), lambda_t)
%         vVspot    Nx1 vector of spot volatility estimates
%         dt        time discretisation
%         vParam    vector of parameters
%         c         double, scaling factor
%
%     Output:
%         out       double, integrand value for given vector of parameters
%
%   author: Evgenii Vladimirov
%   date:   30.03.2019 
%
%% 
    
    mCF = SVhatHJ_ccf2(s, mX, vVspot, dt, vParam, c);
    mExpX = exp(1i.*[c*diff(mX(:,1)), mX(2:end,2)]*s);    
    mH = mExpX - mCF;
    
    iP = size(s,2);
    Int1 = zeros(iP,1);
    parfor j = 1:iP
       HH = mean(exp(1i.*[c*diff(mX(1:end-1,1)), mX(2:end-1,2)]*s(:,j)).*mH(2:end,:),1); 
       Int1(j) = real(HH*(HH'.*WTS));
    end
    
    out = transpose(Int1)*WTS;
end

