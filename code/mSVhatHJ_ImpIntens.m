function [mIntens] = mSVhatHJ_ImpIntens(mOptPrice1, mOptPrice2, mY, mVspot, mK1, mK2, tau, r, mParam)
%   Imply jump intensities for 2SVhatHJ model 
%
%   Used to produce Table R.1, Table R.2, Table R.5, Figure R.1 in the replication pdf 
%       (analogous to Table 3, Table 4, Table 2 and C.1, and Figure 2 in the main text)
%
%     Inputs:
%         mOptPrice1  iNxiKxiT 3d matrix, observed option prices in BSIV
%                                         terms for the 1st stock
%         mOptPrice2  iNxiKxiT 3d matrix, observed option prices in BSIV 
%                                         terms for the 2nd stock
%         mY          iNx2 matrix of current log stock prices
%         mVspot      iNx2 matrix of volatility spot esimates
%         mK1         iNxiKxiT 3d matrix, strike prices, 1st stock
%         mK2         iNxiKxiT 3d matrix, strike prices, 2nd stock
%         tau         double, time-to-maturity, in YEARs
%         r           double, risk-free return, annual
%         mParam      vector of parameters
%
%     Output:
%         mIntens     iNx2 matrix, implied intensities  
%
%   author: Evgenii Vladimirov
%   date:   24.04.2019
%
%%
    [iN,iK,iTau] = size(mOptPrice1);
    mIntens = zeros(iN,2);
    mSolODE1 = mSVhatHJ_ODE_cos(tau, r, mParam);
    mSolODE2 = mSVhatHJ_ODE_cos(tau, r, flip(mParam));
    
    start = mParam(:,4)';
    options = optimset('Display', 'off', 'TolFun', 1e-6, 'TolX', 1e-6);


    % Routine coordinates for COS method:
    a  = -5.0;  %integration lower bound
    b  =  5.0;  %integration upper bound
    c  =  0.0;  %other integration bound
    d  =  b;    %integration bound
    N_cos = 1024;  %Integration bound for COS pricing part
    k = cumsum(ones(N_cos,1))-1; %Integration grid values

    %Routine set-up:
    arg_dma  = k*pi*(d-a)/(b-a);
    arg_cma  = k*pi*(c-a)/(b-a);
    fact_sin = k*pi/(b-a);
    chi = cos(arg_dma)*exp(d) - cos(arg_cma)*exp(c) + fact_sin.*sin(arg_dma)*exp(d) - fact_sin.*sin(arg_cma)*exp(c); 
    chi = chi.*((1 + fact_sin.^2).^(-1));
    phi = (sin(arg_dma)-sin(arg_cma)).*[0; (b-a).*((k(2:end)*pi).^(-1))] + [(d-c);zeros(N_cos-1,1)];
    V_call = 2.0/(b-a)*(chi-phi);
    aux = exp(-1i*(k*a*pi/(b-a))).*V_call;

    parfor t=1:iN
        mSol11 = mSolODE1(:,1,:) + permute(0.5*mVspot(t,1)*(-mSolODE1(:,2,1) + mSolODE1(:,2,1).^2).*tau',[1,3,2]);
        mSol12 = mSolODE2(:,1,:) + permute(0.5*mVspot(t,2)*(-mSolODE2(:,2,1) + mSolODE2(:,2,1).^2).*tau',[1,3,2]);

        [meshTau, meshK1] = meshgrid(tau, mK1(t,:));
        %precalculation for option pricing function
        vCFaux1 = mSol11 + mSolODE1(:,2,:)*mY(t,1);
        aux21 =  1i*k*pi/(b-a).*(log(mK1(t,:)));
        mCFaux1 = bsxfun(@(a,b)(exp(a-b)), vCFaux1, aux21).*aux;
        
        [meshTau, meshK2] = meshgrid(tau, mK2(t,:));
        %precalculation for option pricing function
        vCFaux2 = mSol12 + mSolODE2(:,2,:)*mY(t,2);
        aux22 =  1i*k*pi/(b-a).*(log(mK2(t,:)));
        mCFaux2 = bsxfun(@(a,b)(exp(a-b)), vCFaux2, aux22).*aux;

        func = @(x) reshape([permute(mOptPrice1(t,:,:),[2,3,1]) - GetPrices2(mY(t,1),x(1),x(2), r, mSolODE1(:,3,:),mSolODE1(:,4,:),mCFaux1, meshK1, meshTau);... 
                             permute(mOptPrice2(t,:,:),[2,3,1]) - GetPrices2(mY(t,2),x(2),x(1), r, mSolODE2(:,3,:),mSolODE2(:,4,:),mCFaux2, meshK2, meshTau)],2*iK*iTau,1);
       
        mIntens(t,:) = lsqnonlin(func, start, [0,0], [300, 300], options);
   end    
end

function [out] = GetPrices2(y,lambda1, lambda2, r, mSol_lambda1, mSol_lambda2, mCFaux, meshK, meshTau)
% Calculate option prices using COS method with precalculated parameters
%   
%   Inputs:
%       y               double  log stock price
%       lambda          double  intensity
%       r               double  interest rate
%       mSolODE_lambda  Nx1x4   beta(2) from ODE solution
%       mCFaux          NxKxTau precalculated CF without intensity part
%       meshK
%       meshTau
%
%% Note: the function is modified specifically to back out intensity in SVhatHJ model
%        
%%
    aux3 = mCFaux.*exp(mSol_lambda1*lambda1 + mSol_lambda2*lambda2);
    aux4 = permute(real(sum([0.5*aux3(1,:,:); aux3(2:end,:,:)],1)),[2,3,1]);
    
    out = exp(-r*meshTau).*meshK.*aux4;
    out = calcBSImpVol(1,out,exp(y),meshK,meshTau,r,0);
    
end

