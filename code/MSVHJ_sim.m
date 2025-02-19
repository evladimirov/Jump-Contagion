function [mY_hf, mY, mV, mLambda, count] = MSVHJ_sim(iN, dt, iM, mX0, r, mParam, corr, bPQ)
%   Simulate state data from the 2SVHJ 
%
%       Inputs:
%           iN      integer, number of observations
%           dt      double, time discretisation (dt = iT/iN) between two
%                           ovservations
%           iM      integer, number of draws in between two data points
%           mX0     3xK matrix, initial state vector
%           r       double, risk-free rate
%           mParam  KxP matrix of parameters for each of the K stocks
%           bPQ     boolean, 1 if simulate under P
%                            0 if simulate under Q
%
%       Return:
%           mY_hf       (iN+1)*100x2 matrix of HF stock prices
%           mY          (iN+1)x2 matrix of simulated stock path
%           mV          (iN+1)x2 matrix of simulated volatilities
%           mLambda     (iN+1)x2 matrix of intensities
%
%   author: Evgenii Vladimirov
%   date:   22.03.2019
%   updated: 09.01.2020
%
%% 
    iP = iM*iN;
    dtt = dt/iM;
    
    mdW = sqrt(dtt)*randn(iP,4); %components of Brownian motions
    vB = corr*mdW(:,1)+sqrt(1-corr^2)*mdW(:,3); %brownian component in the secon process
    
    mU = rand(iP,2);
     
    mY = [mX0(1,:); zeros(iP,2)];
    mV = [mX0(2,:); zeros(iP,2)];
    mLambda = [mX0(3,:); zeros(iP,2)];
    
    count1 = 0;
    count2 = 0;
    if bPQ == 1 %simulate states under P
        [sigma_v1, kv1, vbar1, rho1] = deal(mParam(1,1),mParam(1,2), mParam(1,3), mParam(1,4));
        [muj_q1, sigmaj1, kl1, lbar1, delta1, mut_delta1] = deal(mParam(1,5),mParam(1,6),mParam(1,7), mParam(1,8), mParam(1,9), mParam(1,10));
        [muj1, eta1] = deal(mParam(1,11), mParam(1,12));
        
        [sigma_v2, kv2, vbar2, rho2] = deal(mParam(2,1),mParam(2,2), mParam(2,3), mParam(2,4));
        [muj_q2, sigmaj2, kl2, lbar2, delta2, mut_delta2] = deal(mParam(2,5),mParam(2,6),mParam(2,7), mParam(2,8), mParam(2,9), mParam(2,10));
        [muj2, eta2] = deal(mParam(2,11), mParam(2,12));
        
        for i=2:(iP+1)

            mV(i,1) = max(mV(i-1,1) + kv1*(vbar1 - mV(i-1,1))*dtt + sigma_v1*sqrt(mV(i-1,1))*(rho1*mdW(i-1,1) + sqrt(1-rho1^2)*mdW(i-1,2)),0);
            mV(i,2) = max(mV(i-1,2) + kv2*(vbar2 - mV(i-1,2))*dtt + sigma_v2*sqrt(mV(i-1,2))*(rho2*vB(i-1) + sqrt(1-rho2^2)*mdW(i-1,4)),0);
            
            mLambda(i,1) = mLambda(i-1,1) + kl1*(lbar1 - mLambda(i-1,1))*dtt + delta1*(mLambda(i-1,1)*dtt > mU(i-1,1)) + mut_delta1*(mLambda(i-1,2)*dtt > mU(i-1,2));
            mLambda(i,2) = mLambda(i-1,2) + kl2*(lbar2 - mLambda(i-1,2))*dtt + mut_delta2*(mLambda(i-1,1)*dtt > mU(i-1,1)) + delta2*(mLambda(i-1,2)*dtt > mU(i-1,2));

            mY(i,1) = mY(i-1,1) + ((eta1-0.5)*mV(i-1,1)-(exp(muj_q1+0.5*sigmaj1^2)-1)*mLambda(i-1,1))*dtt ...
                + sqrt(mV(i-1,1))*mdW(i-1,1) + (mLambda(i-1,1)*dtt > mU(i-1,1))*(muj1+sigmaj1*randn(1));
            mY(i,2) = mY(i-1,2) + ((eta2-0.5)*mV(i-1,2)-(exp(muj_q2+0.5*sigmaj2^2)-1)*mLambda(i-1,2))*dtt ...
                + sqrt(mV(i-1,2))*vB(i-1) + (mLambda(i-1,2)*dtt > mU(i-1,2))*(muj2+sigmaj2*randn(1));
            
            count1 = count1 + (mLambda(i-1,1)*dtt > mU(i-1,1));
            count2 = count2 + (mLambda(i-1,2)*dtt > mU(i-1,2));
        end
        
    else %simulate under Q 

        [sigma_v1, kv1, vbar1, rho1] = deal(mParam(1,1),mParam(1,2), mParam(1,3), mParam(1,4));
        [muj_q1, sigmaj1, kl1, lbar1, delta1, mut_delta1] = deal(mParam(1,5),mParam(1,6),mParam(1,7), mParam(1,8), mParam(1,9), mParam(1,10));
        [sigma_v2, kv2, vbar2, rho2] = deal(mParam(2,1),mParam(2,2), mParam(2,3), mParam(2,4));
        [muj_q2, sigmaj2, kl2, lbar2, delta2, mut_delta2] = deal(mParam(2,5),mParam(2,6),mParam(2,7), mParam(2,8), mParam(2,9), mParam(2,10));
 
        for i=2:(iP+1)
            mV(i,1) = max(mV(i-1,1) + kv1*(vbar1 - mV(i-1,1))*dtt + sigma_v1*sqrt(mV(i-1,1))*(rho1*mdW(i-1,1) + sqrt(1-rho1^2)*mdW(i-1,2)),0);
            mV(i,2) = max(mV(i-1,2) + kv2*(vbar2 - mV(i-1,2))*dtt + sigma_v2*sqrt(mV(i-1,2))*(rho2*mdW(i-1,3) + sqrt(1-rho2^2)*mdW(i-1,4)),0);
            
            mLambda(i,1) = mLambda(i-1,1) + kl1*(lbar1 - mLambda(i-1,1))*dtt + delta1*(mLambda(i-1,1)*dtt > mU(i-1,1)) + mut_delta1*(mLambda(i-1,2)*dtt > mU(i-1,2));
            mLambda(i,2) = mLambda(i-1,2) + kl2*(lbar2 - mLambda(i-1,2))*dtt + mut_delta2*(mLambda(i-1,1)*dtt > mU(i-1,1)) + delta2*(mLambda(i-1,2)*dtt > mU(i-1,2));

            mY(i,1) = mY(i-1,1) + (r-0.5*mV(i-1,1)-(exp(muj_q1+0.5*sigmaj1^2)-1)*mLambda(i-1,1))*dtt ...
                + sqrt(mV(i-1,1))*mdW(i-1,1) + (mLambda(i-1,1)*dtt > mU(i-1,1))*(muj_q1+sigmaj1*randn(1));
            mY(i,2) = mY(i-1,2) + (r-0.5*mV(i-1,2)-(exp(muj_q2+0.5*sigmaj2^2)-1)*mLambda(i-1,2))*dtt ...
                + sqrt(mV(i-1,2))*mdW(i-1,3) + (mLambda(i-1,2)*dtt > mU(i-1,2))*(muj_q2+sigmaj2*randn(1));
        end
    end
    mY_hf = mY;
    mY = mY(1:iM:iP+1,:);
    mV = mV(1:iM:iP+1,:); 
    mLambda= mLambda(1:iM:iP+1,:);
    count = [count1, count2];
    disp(["Number of jumps:", count])
end

