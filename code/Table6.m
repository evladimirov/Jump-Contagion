function [results] = Table6()
%
%   This function does calculations necessary for Table 6 given the
%   parameter estimates of the bivariate and univariate models for FTSE 100
%   and S&P 500

    
    % Load the estimates reported in the paper
    load('data/estimates.mat', 'mThetaBi', 'mThetaUni', 'r')

    % flip to have the same order in the final table
    mThetaBi = flip(mThetaBi);
    mThetaUni = flip(mThetaUni);

    % Setup necessary parameters for the simulations
    % initial values
    S10 = 100; y10 = log(S10); v1 = 0.007;  
    S20 = 100; y20 = log(S20); v2 = 0.007; 
    
    cor = 0.6;       %contemporaneous correlation between brownians of indices
    
    iN  = 10;        %number of observation
    dt  = 1/365;     %daily observations
    n   = 100;       %number of intraday returns
    iMC = 100000;    %number of simulations


%% simulate 

% fix seedoffset for replicability
seedoffset = iMC;

% Define scenarios: [lambda1_0, lambda2_0]
scenarios = {
    '(a): Base Case',         'base', 'base';
    '(b): Euro Debt Crisis',       5,     5;
    '(c): S&P Shock',         'base',    20; 
    '(d): FTSE Shock',            20, 'base';
    '(e): 2008 Financial Crisis', 15,    20
};

model_list = {'Bivariate - FTSE', 'Univariate - FTSE', 'Bivariate - S&P', 'Univariate - S&P'};

Q = [0.001, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95]; % Quantile levels
results = []; % Store results

disp('Table 6: Descriptive statistics for the conditional log-return distribution (simulated using model parameter estimates, horizon h= 10 days)')
disp(' ')

% Loop over each scenario
for s = 1:size(scenarios, 1)
    scenario_name = scenarios{s,1};
    lambda1_0 = scenarios{s,2};
    lambda2_0 = scenarios{s,3};
    
    
    disp(['                     Scenario ', scenario_name])

    if ischar(lambda1_0)
        lambda1_0 = mThetaBi(1,4);
    end

    if ischar(lambda2_0)
        lambda2_0 = mThetaBi(2,4);
    end

    % Initialize matrices
    simY_bi = zeros(iN+1,2,iMC);
    mCount_bi = zeros(iN+1,2,iMC);
    mX0_bi = [[y10; lambda1_0], [y20; lambda2_0]]; % Adjust initial lambda values

    % Simulate from the bivariate model estimates
    parfor i=1:iMC    
        rng(seedoffset+i, 'twister');
        [~, mY, ~, mCount] = mSVhatHJ_sim(iN, dt, n, mX0_bi, r, mThetaBi, cor, v1, v2);
        simY_bi(:,:,i) = mY;
        mCount_bi(:,:,i) = mCount;
    end

    % Adjust parameters for the univaraite results if necessary
    if lambda1_0 == mThetaBi(1,4)
        lambda1_0 = mThetaUni(1,4);
    end

    if lambda2_0 == mThetaBi(2,4)
        lambda2_0 = mThetaUni(2,4);
    end

    simY_uni = zeros(iN+1,2,iMC);
    mCount_uni = zeros(iN+1,2,iMC);
    mX0_uni = [[y10; lambda1_0], [y20; lambda2_0]]; 

    % Simulate from the univariate model estimates
    parfor i=1:iMC    
        rng(seedoffset+i, 'twister');
        [~, mY, ~, mCount] = mSVhatHJ_sim(iN, dt, n, mX0_uni, r, mThetaUni, cor, v1, v2);
        simY_uni(:,:,i) = mY;
        mCount_uni(:,:,i) = mCount;
    end

    % Compute 10-day log-returns
    mRet10_bi  = squeeze(simY_bi(10,:,:) - [y10, y20]);
    mRet10_uni = squeeze(simY_uni(10,:,:) - [y10, y20]);

    % Compute Quantiles
    Quant10 = 100 * [quantile(mRet10_bi(1,:), Q); quantile(mRet10_uni(1,:), Q);
                     quantile(mRet10_bi(2,:), Q); quantile(mRet10_uni(2,:), Q)];

    % Compute Skewness and Kurtosis
    SK10  = [[skewness(mRet10_bi(1,:)); skewness(mRet10_uni(1,:)); 
              skewness(mRet10_bi(2,:)); skewness(mRet10_uni(2,:))], ...
             [kurtosis(mRet10_bi(1,:)); kurtosis(mRet10_uni(1,:)); 
              kurtosis(mRet10_bi(2,:)); kurtosis(mRet10_uni(2,:))]];

    % Compute Expected Number of Jumps
    numJumps10 = [mean(mCount_bi(10,1,:)); mean(mCount_uni(10,1,:));
                  mean(mCount_bi(10,2,:)); mean(mCount_uni(10,2,:))];

    for i=1:numel(model_list)

        % Store results
        results = [results; {scenario_name, model_list{i},  Quant10(i,:), SK10(i,:), numJumps10(i)}];

        % Print row values with correct format
        fprintf('%-20s | %6.2f | %6.2f | %6.2f | %6.2f | %6.2f | %6.2f | %6.2f | %6.2f | %6.2f | %.4f\n', ...
             model_list{i}, Quant10(i,1), Quant10(i,2), Quant10(i,3), Quant10(i,4), Quant10(i,5), Quant10(i,6), Quant10(i,7), SK10(i,1), SK10(i,2), numJumps10(i));
    end
    disp('-----------------------------------------------------------------------------------------------------------------');

end


end

