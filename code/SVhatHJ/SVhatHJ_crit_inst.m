function [out] = SVhatHJ_crit_inst(theta, PTS, WTS, vY, mOptPrice, vVspot, mK, tau, dt, r, c)
%  Criterion function for C-GMM for SVhatHJ model with check for constraint  
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
%         c         double, scaling factor
%
%     Output:
%         out       double, value of objective function
%
%   author: Evgenii Vladimirov
%   date:   01.04.2019 
%   last version: 08.05.2019
%
%%
    disp(theta)

    % Back-out states    
    [vIntens] = SVhatHJ_ImpIntens(mOptPrice, vY, vVspot, mK, tau, r, theta);

    % calculate moment condition
    out = SVhatHJ_int_inst(PTS, WTS, [vY, vIntens], vVspot, dt, theta, c);
    disp(out)

end

