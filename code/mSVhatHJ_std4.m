function [mA, mB] = mSVhatHJ_std4(theta, PTS, WTS, mY, mOptPrice1, mOptPrice2, mVspot, mK1, mK2, tau, dt, r, c)
%  
%   Std calculation routine for mSVhatHJ
%
%   Used to produce Table R.1 in the replication pdf (analogous to Table 3 in the main text)
%
%     Inputs:
%         theta     vector of parameters
%         PTS       2xP matrix of stacked arguments
%         WTS       1xP vector of weights for each set of argumets
%         vY        (N+1)x1 vector of log stock prices
%         mOptPrice (N+1)xK matrix of option prices
%         vVspot    Nx1 vector of volatility spot estimates based on HF
%         mK        (N+1)xK matrix of stike prices
%         tau       double, fixed time-to-maturity
%         dt        time discretisation
%         r         double, interest rate
%
%     Output:
%         mA        matrix A for the sandwich
%         mB        matrix B for the sandwich
%
%   author: Evgenii Vladimirov
%   date:   16.10.2019 
%
%%
    iK = length(theta);
    iP = size(PTS,2);
    
    %step sizes
    vStep = 1e-3*(abs(theta) + 1e-5);
    vStep = max(vStep', 5e-4);
    mStep = diag(vStep);
    
    %gradient
    disp('calculating the gradient...')
    mfp1 = zeros(iK,iP); mfm1 = zeros(iK,iP);
    mfp2 = zeros(iK,iP); mfm2 = zeros(iK,iP);
    mfp3 = zeros(iK,iP); mfm3 = zeros(iK,iP);
    mfp4 = zeros(iK,iP); mfm4 = zeros(iK,iP);
    for i=1:iK
        [vHp1,vHp2,vHp3,vHp4,~,~,~,~] = mSVhatHJ_mom4(theta + mStep(i,:), PTS, mY, mOptPrice1, mOptPrice2, mVspot, mK1, mK2, tau, dt, r, c);
        [vHm1,vHm2,vHm3,vHm4,~,~,~,~] = mSVhatHJ_mom4(theta - mStep(i,:), PTS, mY, mOptPrice1, mOptPrice2, mVspot, mK1, mK2, tau, dt, r, c);
        mfp1(i,:) = vHp1; mfp2(i,:) = vHp2; mfp3(i,:) = vHp3; mfp4(i,:) = vHp4;
        mfm1(i,:) = vHm1; mfm2(i,:) = vHm2; mfm3(i,:) = vHm3; mfm4(i,:) = vHm4;
    end
    
    mG1 = (mfp1 - mfm1)./(2*vStep); %iK x iP matrix of gradients
    mG2 = (mfp2 - mfm2)./(2*vStep); %iK x iP matrix of gradients
    mG3 = (mfp3 - mfm3)./(2*vStep); %iK x iP matrix of gradients
    mG4 = (mfp4 - mfm4)./(2*vStep); %iK x iP matrix of gradients
    
    disp('calculating A matrix...')
    mA = zeros(iK,iK); %each element is an integral inner product
    for i=1:iP
       mA = mA + real(mG1(:,i)*mG1(:,i)' + mG2(:,i)*mG2(:,i)' + mG3(:,i)*mG3(:,i)' + mG4(:,i)*mG4(:,i)')*WTS(i); 
    end
    
    %% calculate the middle term, matrix B
    disp('calculating moments...')
    [~,~,~,~,mH1,mH2,mH3,mH4] = mSVhatHJ_mom4(theta, PTS, mY, mOptPrice1, mOptPrice2, mVspot, mK1, mK2, tau, dt, r, c);
        
    disp('calculating the covariance operator...')
    mCovOp11 = zeros(iP,iK); mCovOp12 = zeros(iP,iK); mCovOp13 = zeros(iP,iK); mCovOp14 = zeros(iP,iK);
    mCovOp22 = zeros(iP,iK); mCovOp23 = zeros(iP,iK); mCovOp24 = zeros(iP,iK);
    mCovOp33 = zeros(iP,iK); mCovOp34 = zeros(iP,iK);
    mCovOp44 = zeros(iP,iK);
    
    parfor i=1:iP
       %mCovOp(i,:) = (transpose(mH(:,i))*conj(mH))* transpose(mG.*WTS); 
       mCovOp11(i,:) = (transpose(mH1(:,i))*conj(mH1))/size(mH1,1)* transpose(mG1.*WTS); 
       mCovOp12(i,:) = (transpose(mH1(:,i))*conj(mH2))/size(mH2,1)* transpose(mG2.*WTS); 
       mCovOp13(i,:) = (transpose(mH1(:,i))*conj(mH3))/size(mH3,1)* transpose(mG3.*WTS); 
       mCovOp14(i,:) = (transpose(mH1(:,i))*conj(mH4))/size(mH4,1)* transpose(mG4.*WTS); 
       
       mCovOp22(i,:) = (transpose(mH2(:,i))*conj(mH2))/size(mH2,1)* transpose(mG2.*WTS); 
       mCovOp23(i,:) = (transpose(mH2(:,i))*conj(mH3))/size(mH3,1)* transpose(mG3.*WTS); 
       mCovOp24(i,:) = (transpose(mH2(:,i))*conj(mH4))/size(mH4,1)* transpose(mG4.*WTS); 
       
       mCovOp33(i,:) = (transpose(mH3(:,i))*conj(mH3))/size(mH3,1)* transpose(mG3.*WTS); 
       mCovOp34(i,:) = (transpose(mH3(:,i))*conj(mH4))/size(mH4,1)* transpose(mG4.*WTS); 
       
       mCovOp44(i,:) = (transpose(mH4(:,i))*conj(mH4))/size(mH4,1)* transpose(mG4.*WTS); 
    end
    
    disp('calculating B matrix...')
    mB = zeros(iK, iK);
    for i=1:iP
       mB = mB + ...
           (mG1(:,i)*conj(mCovOp11(i,:)) + 2*mG1(:,i)*conj(mCovOp12(i,:)) + 2*mG1(:,i)*conj(mCovOp13(i,:)) + 2*mG1(:,i)*conj(mCovOp14(i,:))...
           + mG2(:,i)*conj(mCovOp22(i,:)) + 2*mG2(:,i)*conj(mCovOp23(i,:)) + 2*mG2(:,i)*conj(mCovOp24(i,:)) ...
           + mG3(:,i)*conj(mCovOp33(i,:)) + 2*mG3(:,i)*conj(mCovOp34(i,:))...
           + mG4(:,i)*conj(mCovOp44(i,:)))*WTS(i); 
    end
    
end


