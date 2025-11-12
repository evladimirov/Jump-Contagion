function [mCF1, mCF2, mCF3, mCF4] = mSVhatHJ_ccf4(s, mLatIntens, mVspot, dt, mParam, c)
%
%     Integrand of the criterion function for the 1st step C-GMM for 2SVhatHJ model
% 
%     Used to produce Table R.1 in the replication pdf (analogous to Table 3 in the main text)
%
%     Inputs:
%         s           4xP matrix of stacked arguments
%         mY          Nx2 matrix of log prices
%         mLatIntens  Nx2 matrix of intensities
%         mVspot      Nx2 matrix of spot volatility estimates
%         dt          time discretisation
%         mParam      matrix of the following parameters:
%                       muj_q, sigmaj, kappa_l, lbar, delta, mut_delta, muj, eta
%         Ft          NxN matrix of instruments
%
%     Output:
%         mCF1        matrix of characteristic functions for state 1
%
%   author: Evgenii Vladimirov
%   date:   30.03.2019 
%%
    [muj_q1, sigmaj1, kl1, lbar1, delta1, mut_delta1, muj1, eta1] = deal(mParam(1,1),mParam(1,2),mParam(1,3), mParam(1,4), mParam(1,5), mParam(1,6), mParam(1,7), mParam(1,8));
    [muj_q2, sigmaj2, kl2, lbar2, delta2, mut_delta2, muj2, eta2] = deal(mParam(2,1),mParam(2,2),mParam(2,3), mParam(2,4), mParam(2,5), mParam(2,6), mParam(2,7), mParam(2,8));
    
    % Solve ODEs 
    % Define CCF system matrices
    K0 = [ 0;
           0;
         kl1*lbar1;
         kl2*lbar2];
    K1 = [ 0,    0,   -c*(exp(muj_q1 + 0.5*sigmaj1^2)-1),   0;
           0,    0,       0,                             -c*(exp(muj_q2 + 0.5*sigmaj2^2)-1);
           0,    0,      -kl1,                            0;
           0,    0,       0,                             -kl2];

    H0 = zeros(4,4);
    H1 = zeros(4,4,4);

    L0 = [0, 0];
    L1 = [0, 0;
          0, 0;
          1, 0;
          0, 1];
     
    JT = @(vBeta)([exp(muj1*vBeta(1)*c +0.5*c^2*sigmaj1^2*vBeta(1)^2 +vBeta(3)*delta1 + vBeta(4)*mut_delta2)-1; 
                   exp(muj2*vBeta(2)*c +0.5*c^2*sigmaj2^2*vBeta(2)^2 +vBeta(3)*mut_delta1 + vBeta(4)*delta2)-1]);

    iP = size(s,2);
    mSolODE_set1 = zeros(5,iP);
    mSolODE_set2 = zeros(5,iP);
    mSolODE_set3 = zeros(5,iP);
    mSolODE_set4 = zeros(5,iP);
    ss1 = [s(1,:);zeros(1,iP);s(2,:);zeros(1,iP)];
    ss2 = [zeros(1,iP);s(1,:);zeros(1,iP);s(2,:)];
    ss3 = [s(1,:);zeros(2,iP);s(2,:)];
    ss4 = [zeros(1,iP);s;zeros(1,iP)];
    parfor j =1:iP
        mSolODE_set1(:,j) = permute(affineODE(ss1(:,j), dt, K0, K1, H0, H1, L0, L1, JT),[1,3,2]);
        mSolODE_set2(:,j) = permute(affineODE(ss2(:,j), dt, K0, K1, H0, H1, L0, L1, JT),[1,3,2]);
        mSolODE_set3(:,j) = permute(affineODE(ss3(:,j), dt, K0, K1, H0, H1, L0, L1, JT),[1,3,2]);
        mSolODE_set4(:,j) = permute(affineODE(ss4(:,j), dt, K0, K1, H0, H1, L0, L1, JT),[1,3,2]);
    end
    
    % Calculate four Characteristic Functions
    mCF1 = exp((mVspot(1:end-1,1)*((eta1-0.5)*1i*ss1(1,:)*c - 0.5*ss1(1,:).^2*c^2))*dt...
        + mSolODE_set1(1,:) + mLatIntens(1:end-1,1)*mSolODE_set1(4,:) + mLatIntens(1:end-1,2)*mSolODE_set1(5,:));
    mCF2 = exp((mVspot(1:end-1,2)*((eta2-0.5)*1i*ss2(2,:)*c - 0.5*ss2(2,:).^2*c^2))*dt...
        + mSolODE_set2(1,:) + mLatIntens(1:end-1,1)*mSolODE_set2(4,:) + mLatIntens(1:end-1,2)*mSolODE_set2(5,:));
    
    mCF3 = exp((mVspot(1:end-1,1)*((eta1-0.5)*1i*ss1(1,:)*c - 0.5*ss1(1,:).^2*c^2))*dt...
        + mSolODE_set3(1,:) + mLatIntens(1:end-1,1)*mSolODE_set3(4,:) + mLatIntens(1:end-1,2)*mSolODE_set3(5,:));
    mCF4 = exp((mVspot(1:end-1,2)*((eta2-0.5)*1i*ss2(2,:)*c - 0.5*ss2(2,:).^2*c^2))*dt...
        + mSolODE_set4(1,:) + mLatIntens(1:end-1,1)*mSolODE_set4(4,:) + mLatIntens(1:end-1,2)*mSolODE_set4(5,:));
    
end
