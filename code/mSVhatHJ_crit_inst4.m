function [out] = mSVhatHJ_crit_inst4(theta, PTS, WTS, mY, mOptPrice1, mOptPrice2, mVspot, mK1, mK2, tau, dt, r, c)
%  Criterion function for C-GMM for 2SVhatHJ model with constraint  
% 
%     Used to produce Table R.1 and Table R.5 in the replication pdf 
%       (analogous to Table 3 and Table 2 and C.1 in the main text)
%
%     Inputs:
%         theta         vector of parameters
%         PTS           2xP matrix of stacked arguments
%         WTS           1xP vector of weights for each set of argumets
%         mY            (N+1)x2 vector of log stock prices
%         mOptPrice1    (N+1)xK matrix of option prices on stock 1
%         mOptPrice2    (N+1)xK matrix of option prices on stock 2
%         mVspot        Nx2 vector of volatility spot estimates based on HF
%         mK1           (N+1)xK matrix of stike prices for stock 1
%         mK2           (N+1)xK matrix of stike prices for stock 2
%         tau           double, fixed time-to-maturity
%         dt            time discretisation
%         r             double, interest rate
%         c             double, scaling factor
%
%     Output:
%         out           double, value of objective function
%
%   author: Evgenii Vladimirov
%   date:   25.04.2019 
%   last version: 08.05.2019
%
%%
    mTheta = [theta(1:end/2); theta(end/2+1:end)];
    disp(mTheta)
        
    % Back-out intensities
    mIntens = mSVhatHJ_ImpIntens(mOptPrice1(2:end,:,:), mOptPrice2(2:end,:,:), mY(2:end,:), mVspot, mK1(2:end,:), mK2(2:end,:), tau, r, mTheta);

    % calculate moment condition
    [out1, out2, out3, out4] = mSVhatHJ_int_inst4(PTS, WTS, mY(2:end,:), mIntens, mVspot, dt, mTheta, c);
    out = out1 + out2 + out3 + out4;
    disp([out1, out2, out3, out4, out])

end

