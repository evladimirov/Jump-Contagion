function [vH1, vH2, vH3, vH4, mH1, mH2, mH3, mH4] = mSVhatHJ_mom4(theta, PTS, mY, mOptPrice1, mOptPrice2, mVspot, mK1, mK2, tau, dt, r, c)
%
%   Calculate values of the moment function for a given value of parameters
%
%   Used to produce Table R.1 in the replication pdf (analogous to Table 3 in the main text)
%
%     Inputs:
%         theta     vector of parameters
%         PTS       2xP matrix of stacked arguments
%         vY        (N+1)x1 vector of log stock prices
%         mOptPrice (N+1)xK matrix of option prices
%         vVspot    Nx1 vector of volatility spot estimates based on HF
%         mK        (N+1)xK matrix of stike prices
%         tau       double, fixed time-to-maturity
%         dt        time discretisation
%         r         double, interest rate
%
%     Output:
%         vH        1xP vector of moment values
%         mH        NxP marix of moment values, vH=sum(mH)
%
%   author: Evgenii Vladimirov
%   date:   29.06.2020 
%
%%
    mTheta = [theta(1:8); theta(9:end)];
    % Back-out states    
    mIntens = mSVhatHJ_ImpIntens(mOptPrice1, mOptPrice2, mY, mVspot, mK1, mK2, tau, r, mTheta);
    
    % calculate moment condition
    [mCF1, mCF2, mCF3, mCF4] = mSVhatHJ_ccf4(PTS(3:4,:), mIntens(2:end,:), mVspot(2:end,:), dt, mTheta, c);
    mExpX1 = exp(1i.*[c*diff(mY(2:end,1)), mIntens(3:end,1)]*PTS(3:4,:)); 
    mExpX2 = exp(1i.*[c*diff(mY(2:end,2)), mIntens(3:end,2)]*PTS(3:4,:)); 
    mExpX3 = exp(1i.*[c*diff(mY(2:end,1)), mIntens(3:end,2)]*PTS(3:4,:)); 
    mExpX4 = exp(1i.*[c*diff(mY(2:end,2)), mIntens(3:end,1)]*PTS(3:4,:)); 
    
    mInstr1 = exp(1i.*[c*diff(mY(1:end-1,1)), mIntens(2:end-1,1)]*PTS(1:2,:));
    mInstr2 = exp(1i.*[c*diff(mY(1:end-1,2)), mIntens(2:end-1,2)]*PTS(1:2,:));
    mInstr3 = exp(1i.*[c*diff(mY(1:end-1,1)), mIntens(2:end-1,2)]*PTS(1:2,:));
    mInstr4 = exp(1i.*[c*diff(mY(1:end-1,2)), mIntens(2:end-1,1)]*PTS(1:2,:));
    
    mH1 = mInstr1.*(mExpX1 - mCF1);
    vH1 = sum(mH1)/(size(mH1,1));
    
    mH2 = mInstr2.*(mExpX2 - mCF2);
    vH2 = sum(mH2)/(size(mH2,1));
    
    mH3 = mInstr3.*(mExpX3 - mCF3);
    vH3 = sum(mH3)/(size(mH3,1));
    
    mH4 = mInstr4.*(mExpX4 - mCF4);
    vH4 = sum(mH4)/(size(mH4,1));
    
end

