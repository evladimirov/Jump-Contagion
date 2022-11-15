function [c, ceq] = mSVhatHJ_fmin_constr(x)
%   2SVhatHJ constraint check for fmincon
%
%   Input:
%       mParam  matrix of parameters
%
%   Output:
%       out     boolean, 1 if constraints are satisfied
%                        0 if not
%
%   author: Evgenii Vladimirov
%   date:   08.05.2019 
%
%%
%   k1 = x(3);  a11 = x(5); a12 = x(6);
%   k2 = x(11); a22= x(13); a21 = x(14);
%
%   k1 = mParam(1,3); a11 = mParam(1,5); a12 = mParam(1,6);
%   k2 = mParam(2,3); a22 = mParam(2,5); a21 = mParam(2,6);
%
%   eigen1 = 0.5*(a22/k2 + a11/k1 + sqrt((a22/k2 - a11/k1)^2 + 4*a21*a12/(k1*k2)));
%   eigen2 = 0.5*(a22/k2 + a11/k1 - sqrt((a22/k2 - a11/k1)^2 + 4*a21*a12/(k1*k2)));
%
    c(1) = abs(0.5*(x(13)/x(11) + x(5)/x(3) + sqrt((x(13)/x(11) - x(5)/x(3))^2 + 4*x(14)*x(6)/(x(3)*x(11)))))-1+eps;
    c(2) = abs(0.5*(x(13)/x(11) + x(5)/x(3) - sqrt((x(13)/x(11) - x(5)/x(3))^2 + 4*x(14)*x(6)/(x(3)*x(11)))))-1+eps;

    
   ceq = [];
   %out = [abs(eigen1)-1; abs(eigen2)-1];
end