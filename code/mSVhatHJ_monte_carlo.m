function [MCTableRes] = mSVhatHJ_monte_carlo()
    %
    %   Monte-Carlo simulation for bivariate 2SVhatHJ estimation procedure
    %
    %   Used to produce Table R.5 in the replication pdf (Table 2 and C.1 in the main text)
    %
    
    %% DGP parameters
    r = 0.0;
    %        [muj_q, sigmaj, kl, lbar, delta, mut_delta, muj,  eta]
    vParam1 = [-0.13, 0.03, 6.0,   1,    3.0,   1.0,    -0.04, 2.0];
    vParam2 = [-0.13, 0.03, 5.0,   1,    2.0,   3.0,    -0.04, 2.0];
    
    cor = 0.6; %contemporaneous correlation between brownians of indices
    
    % volatility parameters used to simulate SV part
    %           [sigma_v, kv, vbar, rho]
    vParamVol1 = [0.22, 4.8, 0.015, -0.6];
    vParamVol2 = [0.22, 4.8, 0.015, -0.6];
    
    % initial values
    S10 = 100; y10 = log(S10); v10 = 0.01; lambda10 = vParam1(4); X10 = [y10; v10; lambda10];
    S20 = 80;  y20 = log(S20); v20 = 0.01; lambda20 = vParam2(4); X20 = [y20; v20; lambda20];
    
    mX0 = [X10, X20];
    mParam_all = [vParamVol1, vParam1; vParamVol2, vParam2];
    mParam = [vParam1; vParam2];
    vParam = [vParam1, vParam2]; 
    
    
    iN = 1500;      %number of observation
    dt = 1/365;     %daily observations
    n = 100;        %number of intraday returns
    
    % option prices characteristics
    vM = [0.8:0.05:1.15];      % moneyness level
    vTau = [0.1:0.1:0.1]';     % time-to-maturity
    
    %% estimation parameters setting
    
    %vParam1 =   [ mujq, sigmaj, kl,  lbar,    a11,  a12,    muj, eta];
    ub = repmat([ -0.05,  0.09, 40,     5,     30,   30,  -0.01,   8], 1,2);
    lb = repmat([ -0.25, 0.005,  1,   0.05, 0.005,    0,   -0.1,  -5], 1,2);
    
    [ Int, WTS, PTS, INTCLS ] = fwtpts( 2, 11, 'Norm');
    
    idK = [1:7]; idTau = [1:1];
    
    optSearch = optimset('Display','iter', 'MaxIter',500, 'PlotFcns', @optimplotfval);
    
    constr = @(theta)(mSVhatHJ_fmin_constr(theta));
    
    
    %% MC simulation
    iMC = 100;
    mTheta_MC = zeros(iMC, 16);
    vSeeds = load('data/seeds', 'vSeeds');
    vSeeds = vSeeds.vSeeds;
    count = [0,0];
    
    for i=1:iMC
        rng(vSeeds(i));
        % generate starting value +- 5% of the true values
        vStart = 0.2*vParam.*rand(1,16) + vParam - 0.1*vParam;
    
        %simulate data
        while or(count(1)<15, count(2)<15)
            [~, mY, mV, mLambda, count] = MSVHJ_sim(iN, dt, n, mX0, r, mParam_all, cor, 1);
    
            [mOptPriceIV1, mOptPriceIV2] = MSVHJ_sim_opt(iN, mY, mV, mLambda, r, mParam_all, vTau, vM, 1);
            mK1 = vM.*exp(mY(:,1)); mK2 = vM.*exp(mY(:,2));
        end
        
        % estimation
        crit4 = @(theta)(mSVhatHJ_crit_inst4(theta,...
            PTS, WTS', mY,mOptPriceIV1(:,idK,idTau), mOptPriceIV2(:,idK,idTau), mV(2:end,:), ...
            mK1(:,idK), mK2(:,idK), vTau(idTau),dt,r,50));
    
        [vTheta, ~] = fminsearchcon(crit4,vStart,lb,ub,[],[],constr,optSearch);
        
        mTheta_MC(i,:) = vTheta;
        
        count = [0,0];
    end
    
    %%

    % Combine statistical summaries into one matrix
    MCres = [
        vParam;                      
        mean(mTheta_MC);               
        sqrt(mean((mTheta_MC - mean(mTheta_MC)).^2));   
        quantile(mTheta_MC, 0.25);  
        quantile(mTheta_MC, 0.5);  
        quantile(mTheta_MC, 0.75)   
    ];

    % Convert to a table with row names and variable names
    MCTableRes = array2table(round(MCres,3), ...
        "RowNames", {'true', 'mean', 'std', '25%', '50%', '75%'}, ...
        "VariableNames", {'mu1_q1', 'sigma1', 'kappa1', 'l_bar1', 'd11', 'd12', 'mu1', 'eta1',...
                          'mu2_q2', 'sigma2', 'kappa2', 'l_bar2', 'd22', 'd21', 'mu2', 'eta2'})

end

