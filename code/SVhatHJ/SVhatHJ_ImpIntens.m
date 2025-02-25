function [vIntens] = SVhatHJ_ImpIntens(mOptPrice, vY, vVspot, mK, tau, r, vParam)
%   Imply latent intensities for SVhatHJ model using option panel and spot volatility
%   
%     Inputs:
%         mOptPrice   iNxiK matrix, observed option prices in BSIV terms
%         vY          iNx1 vector of current log stock prices
%         vVspot      iNx1 vector of spot volatility estimates
%         mK          iNxiK matrix, strike prices
%         tau         iTau x1 vector, time-to-maturity, in YEARs
%         r           double, risk-free return, annual
%         vParam      vector of parameters
%
%     Output:
%         vIntens     iNx1 vector, implied intensities  
%
%   author: Evgenii Vladimirov
%   date:   29.03.2019
%   last version: 05.05.2019
%
%%
    iN = size(mOptPrice,1);
    vIntens = zeros(iN,1);
    mSolODE = SVhatHJ_ODE_cos(tau, r, vParam);
         
    % Routine coordinates for COS method:
    a  = -5.0;  %integration lower bound
    b  =  5.0;  %integration upper bound
    c  =  0.0;  %other integration bound
    d  =  b;    %integration bound
    N_cos 	= 1024;  %Integration bound for COS pricing part
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

    lbar = vParam(4);
    options = optimset('Display', 'off', 'TolFun', 1e-8, 'TolX', 1e-8, 'FinDiffType', 'central');
    
    parfor t=1:iN
        % note this works because beta(1)=u for any t
        mSol1 = mSolODE(:,1,:) + permute(0.5*vVspot(t)*(-mSolODE(:,2,1) + mSolODE(:,2,1).^2).*tau',[1,3,2]);

        [meshTau, meshK] = meshgrid(tau, mK(t,:));
        %precalculations for option pricing function
        vCFaux = mSol1 + mSolODE(:,2,:)*vY(t);
        aux2 =  1i*k*pi/(b-a).*(log(mK(t,:)));
        mCFaux = bsxfun(@(a,b)(exp(a-b)), vCFaux, aux2).*aux;

        func = @(x) sum(permute(mOptPrice(t,:,:),[2,3,1]) - GetPrices2(vY(t),x, r, mSolODE(:,3,:),mCFaux,meshK,meshTau), 2);     
        vIntens(t) = lsqnonlin(func, lbar, 0 , 200, options);
    end    
end


function [out] = GetPrices2(y,lambda, r, mSolODE_lambda, mCFaux, meshK, meshTau)
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
%        for COS method in general see Opt_price_cos.m
%%
    aux3 = mCFaux.*exp(mSolODE_lambda*lambda);
    aux4 = permute(real(sum([0.5*aux3(1,:,:); aux3(2:end,:,:)],1)),[2,3,1]);
    
    out = exp(-r*meshTau).*meshK.*aux4;
    out = calcBSImpVol(1,out,exp(y),meshK,meshTau,r,0);
end



