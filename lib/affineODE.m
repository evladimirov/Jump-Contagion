function [sol] = affineODE(u, T, K0, K1, H0, H1, L0, L1, JT)
% ODE solver for affine model following Duffie, Pan, Singelton (2002)
% 
%     Inputs:
%         u       kx1 vector of arguments
%         T       iTx1 vector of "time-to-maturity", in YEARs
%         K0      kx1 vector
%         K1      kxk matrix
%         H0      kxk matrix
%         H1      kxkxk matrix
%         L0      1xJ vector,
%         L1      kxJ matrix
%         JT      function, "jump transform", which returns Jx1 vector
%     Return:
%         sol     1d vector, solution of the ODE
%
%   author: Evgenii Vladimirov
%   date:   11.03.2019
%   last version: 26.04.2019
%
%% Note:    the relative error of my_ode45 is fixed at 1e-07
%
%%
    func = @(t,vX) diff_eq(t,vX,K0,K1,H0,H1,L0,L1,JT);
    u0 = [0; 1i*u];
    [~,mX] = my_ode45(func,[0;T],u0);
    sol = mX(2:end,:).';
end

%% Creates the RHS of ODE
function vRHS = diff_eq(t,vX,K0,K1,H0,H1,L0,L1,JT)
%
%   Creates the RHS of ODE equation following DPS (2001)
%
%%
    vBeta = vX(2:end);
    aux = zeros(length(vBeta),1);
    for i=1:length(vBeta)
        aux(i,1)=0.5*(vBeta.')*H1(:,:,i)*vBeta;
    end
    
    vRHS = [K0'*vBeta + 0.5*(vBeta.')*H0*vBeta + L0*JT(vBeta);
            K1'*vBeta + aux + L1*JT(vBeta)];
end

