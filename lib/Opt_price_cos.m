function [out] = Opt_price_cos(vX, K, T, r, vParam, bSol, SolODE, bBSIV)
%   Valuation of European call option using COS method of Fang Oosterlee
%     (2008) and numerically solving ODE for CF
% 
%     Inputs:
%         vX        column-vector of states (i.e. for SVJ vX=(log(S0),v0)')
%         K         vector, strike price
%         T         vector, time-to-maturity, in YEARs
%         r         double, annual interest rate
%         vParam    vector of the parameters 
%         bSol      boolean, indicator for pre-calculation of ODE solution
%                   0 = no pre-calculations are provided,                           
%                   1 = pre-calculations are provided
%         SolODE    if bSol=0 -> SolODE has to be passed as a function
%                   if bSol=1 -> SolODE is treated as a matrix of solutions
%         bBSIV     boolean, indicator of the output price
%                   if 0 -> option price in currency value
%                   if 1 -> option price in BS implied volatility terms
%     Return:
%         out       double, option call price
%
%   date:   29.03.2019
%   last version: 26.04.2019
%
%% Note:    function does not allow time-series of state-vector
%%
    % Routine coordinates:
    a  = -5.0; %integration lower bound
    b  =  5.0; %integration upper bound
    c  =  0.0; %other integration bound
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

    if bSol==0
          mSolODE = SolODE(T, r, vParam);
    else
          mSolODE = SolODE;
    end
    
    out = zeros(length(K), length(T));
    for tt=1:length(T)
        for kk=1:length(K)
            vCF = exp(mSolODE(:,1,tt) + mSolODE(:,2:end,tt)*vX - 1i*k*pi/(b-a)*(log(K(kk))));
            out(kk,tt) = exp(-r*T(tt))*K(kk)*real(sum([0.5;ones(N_cos-1,1)].*vCF.*exp(-1i*(k*a*pi/(b-a))).*V_call));
        end
    end 
    
    if bBSIV == 1
%         for tt=1:length(T)
%             for kk=1:length(K)
%                 out(kk,tt) = BS_iv(out(kk,tt), exp(vX(1)), K(kk), T(tt), r); %BS implied volatility
%             end
%         end
        [meshK,meshTau] = meshgrid(K, T);
        out = calcBSImpVol(1,out,exp(vX(1)),meshK',meshTau',r,0);
        %blsimpv(exp(vX(1)), meshK, 0, meshTau, out)
    end

end
