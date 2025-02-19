function Figure3()
%
%   This function does calculations necessary for Figure 3 given the
%   parameter estimates of the bivariate and univariate models for FTSE 100
%   and S&P 500
    
        
    % Load the estimates reported in the paper
    load('data/estimates.mat', 'mThetaBi', 'mThetaUni', 'r')

    % Setup necessary parameters for the simulations
    % initial values
    S10 = 100; y10 = log(S10); v1 = 0.007;  
    S20 = 100; y20 = log(S20); v2 = 0.007; 
    
    cor = 0.6;       %contemporaneous correlation between brownians of indices
    
    iN  = 10;        %number of observation
    dt  = 1/365;     %daily observations
    n   = 100;       %number of intraday returns
    iMC = 100000;    %number of simulations

    % fix seedoffset for replicability
    seedoffset = iMC;

    % Initialize matrices
    simY_bi = zeros(iN+1,2,iMC);
    mX0_bi = [[y10;  20], [y20; mThetaBi(2,4)]]; % Adjust initial lambda values

    % Simulate from the bivariate model estimates
    parfor i=1:iMC    
        rng(seedoffset+i, 'twister');
        [~, mY, ~, ~] = mSVhatHJ_sim(iN, dt, n, mX0_bi, r, mThetaBi, cor, v1, v2);
        simY_bi(:,:,i) = mY;
    end


    simY_uni = zeros(iN+1,2,iMC);
    mX0_uni = [[y10;  20], [y20; mThetaUni(2,4)]]; 

    % Simulate from the univariate model estimates
    parfor i=1:iMC    
        rng(seedoffset+i, 'twister');
        [~, mY, ~, ~] = mSVhatHJ_sim(iN, dt, n, mX0_uni, r, mThetaUni, cor, v1, v2);
        simY_uni(:,:,i) = mY;
    end

    % Compute 10-day log-returns
    mRet10_bi  = squeeze(simY_bi(10,:,:) - [y10, y20]);
    mRet10_uni = squeeze(simY_uni(10,:,:) - [y10, y20]);

    %% plot Figure 3

    y1 = -0.25:0.005:0.2; y2 = -0.4:0.005:0.2;
    pts = zeros(length(y1)*2, 2);
    k=1;
    for i = 1:length(y1)
        for j =1:length(y2)
            pts(k,:) = [y1(i), y2(j)];
            k= k+1;
        end
    end
    
    % Make a plot for the bivariate results
    [f,xi] = ksdensity(mRet10_bi', pts, 'Bandwidth',0.01);  
    
    x1 = unique(xi(:,1)); x2 = unique(xi(:,2));
    Cx = zeros(length(x1), length(x2));
    for i=1:length(x1)
        for j =1:length(x2)
            Cx(i,j) = f(and(xi(:,1)==x1(i), xi(:,2)==x2(j)));
        end
    end
    

    figure
    scatter(mRet10_bi(1,:), mRet10_bi(2,:), '.', 'MarkerEdgeColor', [0.75 0.75 0.75]); hold on
    contour(x1,x2,Cx'); hold on
    xlabel('S&P'); ylabel('FTSE');
    ylim([-0.18, 0.1]);
    xlim([-0.3, 0.25]);

       % add extra contours
    Ck = contourc(x1,x2,Cx',[15, 15]);
    plot(Ck(1,2:end), Ck(2,2:end),'k-','LineWidth',1)
    
    Ck = contourc(x1,x2,Cx',[5, 5]);
    plot(Ck(1,2:end), Ck(2,2:end),'k-','LineWidth',1)
    
    Ck = contourc(x1,x2,Cx',[2, 2]);
    plot(Ck(1,2:end), Ck(2,2:end),'k-','LineWidth',1)
    
    Ck = contourc(x1,x2,Cx',[0.5, 0.5]);
    plot(Ck(1,2:end), Ck(2,2:end),'k-','LineWidth',1)
    
    Ck = contourc(x1,x2,Cx',[0.15, 0.15]);
    plot(Ck(1,2:end), Ck(2,2:end),'k-','LineWidth',1)
    
    %Ck = contourc(x1,x2,Cx',[0.1, 0.1]);
    %plot(Ck(1,2:end), Ck(2,2:end),'k-','LineWidth',1)

    title('(a) Bivariate estimates')
    hold off;


    % Make a plot for the univariate results
    [f,xi] = ksdensity(mRet10_uni', pts, 'Bandwidth',0.01);  

    x1 = unique(xi(:,1)); x2 = unique(xi(:,2));
    Cx = zeros(length(x1), length(x2));
    for i=1:length(x1)
        for j =1:length(x2)
            Cx(i,j) = f(and(xi(:,1)==x1(i), xi(:,2)==x2(j)));
        end
    end
    
    figure
    scatter(mRet10_uni(1,:), mRet10_uni(2,:), '.', 'MarkerEdgeColor', [0.75 0.75 0.75]); hold on
    contour(x1,x2,Cx'); hold on
    xlabel('S&P'); ylabel('FTSE');
    ylim([-0.18, 0.1]);
    xlim([-0.3, 0.25]);

       % add extra contours
    Ck = contourc(x1,x2,Cx',[15, 15]);
    plot(Ck(1,2:end), Ck(2,2:end),'k-','LineWidth',1)
    
    Ck = contourc(x1,x2,Cx',[5, 5]);
    plot(Ck(1,2:end), Ck(2,2:end),'k-','LineWidth',1)
    
    Ck = contourc(x1,x2,Cx',[2, 2]);
    plot(Ck(1,2:end), Ck(2,2:end),'k-','LineWidth',1)
    
    Ck = contourc(x1,x2,Cx',[0.5, 0.5]);
    plot(Ck(1,2:end), Ck(2,2:end),'k-','LineWidth',1)
    
    Ck = contourc(x1,x2,Cx',[0.15, 0.15]);
    plot(Ck(1,2:end), Ck(2,2:end),'k-','LineWidth',1)
    
    title('(b) Univariate estimates')
    hold off;


end

