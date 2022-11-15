function varargout = my_ode45(odeFcn,tspan,y0)
%MY_ODE45  Solve non-stiff differential equations, medium order method.
%   [TOUT,YOUT] = ODE45(ODEFUN,TSPAN,Y0) with TSPAN = [T0 TFINAL] integrates 
%   the system of differential equations y' = f(t,y) from time T0 to TFINAL 
%   with initial conditions Y0. ODEFUN is a function handle. For a scalar T
%   and a vector Y, ODEFUN(T,Y) must return a column vector corresponding 
%   to f(t,y). Each row in the solution array YOUT corresponds to a time 
%   returned in the column vector TOUT.  To obtain solutions at specific 
%   times T0,T1,...,TFINAL (all increasing or all decreasing), use TSPAN = 
%   [T0 T1 ... TFINAL].     

%   A. Lalu, 10.01.13
%   THIS FUNCTION IS A CAUTIOUS SIMPLIFICATION OF THE BUILT-IN ODE SOLVER 
%   AVAILABLE IN THE MATLAB ODE SUITE (R2012b). ALL COPYRIGHTS ARE ATTRIBUTED
%   TO THE ORIGINAL COPYRIGHT HOLDERS AS DETAILED HEREFORTH. 

%   ODE45 is an implementation of the explicit Runge-Kutta (4,5) pair of
%   Dormand and Prince called variously RK5(4)7FM, DOPRI5, DP(4,5) and DP54.
%   It uses a "free" interpolant of order 4 communicated privately by
%   Dormand and Prince.  Local extrapolation is done.

%   Details are to be found in The MATLAB ODE Suite, L. F. Shampine and
%   M. W. Reichelt, SIAM Journal on Scientific Computing, 18-1, 1997.

%   Mark W. Reichelt and Lawrence F. Shampine, 6-14-94
%   ORIGINAL Copyright: 1984-2011 The MathWorks, Inc.
%   $Revision: 5.74.4.13 $  $Date: 2011/04/16 06:38:58 $

% Stats
nsteps  = 0;
nfailed = 0;
nfevals = 0; 

% Set solver arguments
htspan = abs(tspan(2)-tspan(1)); 
ntspan = length(tspan); 
t0 = tspan(1); 
next = 2;
tfinal = tspan(end); 
neq = length(y0); 
tdir = sign(tfinal-t0);
f0 = feval(odeFcn,t0,y0); 
dataType = superiorfloat(t0,y0,f0); 
% rtol = 1e-3;
% atol = 1e-6; 
rtol = 1e-9;
atol = 1e-9; 
threshold = atol/rtol;
hmax = 0.1*(tfinal-t0); 

% Handle the output
if ntspan > 2
  outputAt = 'RequestedPoints';         % output only at tspan points
else
  outputAt = 'SolverSteps';             % computed points, no refinement
end

t = t0;
y = y0;

% Allocate memory if we're generating output.
tout = []; yout = [];
nout = 1;
tout(nout) = t;
yout(:,nout) = y;  

% Initialize method parameters.
pow = 1/5;
A = [1/5, 3/10, 4/5, 8/9, 1, 1];
B = [
    1/5         3/40    44/45   19372/6561      9017/3168       35/384
    0           9/40    -56/15  -25360/2187     -355/33         0
    0           0       32/9    64448/6561      46732/5247      500/1113
    0           0       0       -212/729        49/176          125/192
    0           0       0       0               -5103/18656     -2187/6784
    0           0       0       0               0               11/84
    0           0       0       0               0               0
    ];
E = [71/57600; 0; -71/16695; 71/1920; -17253/339200; 22/525; -1/40];
f = zeros(neq,7,dataType);
hmin = 16*eps(t);
  % Compute an initial step size h using y'(t).
  absh = min(hmax, htspan);
  rh = norm(f0 ./ max(abs(y),threshold),inf) / (0.8 * rtol^pow);
  if absh * rh > 1
    absh = 1 / rh;
  end
  absh = max(absh, hmin);
f(:,1) = f0;


% THE MAIN LOOP
done = false;
while ~done
  
  % By default, hmin is a small number such that t+hmin is only slightly
  % different than t.  It might be 0 if t is 0.
  hmin = 16*eps(t);
  absh = min(hmax, max(hmin, absh));    % couldn't limit absh until new hmin
  h = tdir * absh;
  
  % Stretch the step if within 10% of tfinal-t.
  if 1.1*absh >= abs(tfinal - t)
    h = tfinal - t;
    absh = abs(h);
    done = true;
  end
  
  % LOOP FOR ADVANCING ONE STEP.
  nofailed = true;                      % no failed attempts
  while true
    hA = h * A;
    hB = h * B;
    f(:,2) = feval(odeFcn,t+hA(1),y+f*hB(:,1));
    f(:,3) = feval(odeFcn,t+hA(2),y+f*hB(:,2));
    f(:,4) = feval(odeFcn,t+hA(3),y+f*hB(:,3));
    f(:,5) = feval(odeFcn,t+hA(4),y+f*hB(:,4));
    f(:,6) = feval(odeFcn,t+hA(5),y+f*hB(:,5));

    tnew = t + hA(6);
    if done
      tnew = tfinal;   % Hit end point exactly.
    end
    h = tnew - t;      % Purify h.     
    
    ynew = y + f*hB(:,6);
    f(:,7) = feval(odeFcn,tnew,ynew);         
    
    % Estimate the error.
    err = absh * norm((f * E) ./ max(max(abs(y),abs(ynew)),threshold),inf);  
    % Accept the solution only if the weighted error is no more than the
    % tolerance rtol.  Estimate an h that will yield an error of rtol on
    % the next step or the next try at taking this step, as the case may be,
    % and use 0.8 of this value to avoid failures.
    if err > rtol                       % Failed step
      nfailed = nfailed + 1;            
      if absh <= hmin
        warning(message('MATLAB:ode45:IntegrationTolNotMet', sprintf( '%e', t ), sprintf( '%e', hmin )));
        varargout{1}=tout(1:nout).';
        varargout{2}=yout(:,1:nout).';
        return;
      end
      
      if nofailed
        nofailed = false;



          absh = max(hmin, absh * max(0.1, 0.8*(rtol/err)^pow));

      else
        absh = max(hmin, 0.5 * absh);
      end
      h = tdir * absh;
      done = false;
      
    else                                % Successful step
      break;
    end
  end
  nsteps = nsteps + 1;                  
      

    switch outputAt
     case 'SolverSteps'        % computed points, no refinement
      nout_new=0;
      if tnew == tfinal
         tout = [tout, tnew];
         yout = [yout, ynew];
      end
     case 'RequestedPoints'    % output only at tspan points
      nout_new =  0;
      tout_new = [];
      yout_new = [];
      while next <= ntspan  
        if tdir * (tnew - tspan(next)) < 0
           break;
        end
        nout_new = nout_new + 1;              
        tout_new = [tout_new, tspan(next)];
        if tspan(next) == tnew
          yout_new = [yout_new, ynew];            
        else  
          yout_new = [yout_new, my_ntrp45(tspan(next),t,y,[],[],h,f)];
        end  
        next = next + 1;
      end
      tout = [tout, tout_new];
      yout = [yout, yout_new];   
    end
    nout=nout+nout_new;
  if done
    break
  end

  % If there were no failures compute a new h.
  if nofailed
    % Note that absh may shrink by 0.8, and that err may be 0.
    temp = 1.25*(err/rtol)^pow;
    if temp > 0.2
      absh = absh / temp;
    else
      absh = 5.0*absh;
    end
  end
  
  % Advance the integration one step.
  t = tnew;
  y = ynew;
  f(:,1) = f(:,7);  % Already have f(tnew,ynew)
  
end

varargout{1} = tout.'; 
varargout{2} = yout.';

end



function [yinterp] = my_ntrp45(tinterp,t,y,~,~,h,f)
%NTRP45  Interpolation helper function for ODE45.
%   YINTERP = NTRP45(TINTERP,T,Y,TNEW,YNEW,H,F,IDX) uses data computed in ODE45
%   to approximate the solution at time TINTERP.  TINTERP may be a scalar 
%   or a row vector. 
%   The arguments TNEW and YNEW do not affect the computations. They are 
%   required for consistency of syntax with other interpolation functions. 
%   Any values entered for TNEW and YNEW are ignored.
%    
%   [YINTERP,YPINTERP] = NTRP45(TINTERP,T,Y,TNEW,YNEW,H,F,IDX) returns also the
%   derivative of the polynomial approximating the solution. 
%
%   IDX has indices of solution components that must be non-negative. Negative 
%   YINTERP(IDX) are replaced with zeros and the derivative YPINTERP(IDX) is 
%   set to zero.
%   
%   See also ODE45, DEVAL.

%   Mark W. Reichelt and Lawrence F. Shampine, 6-13-94
%   Copyright 1984-2009 The MathWorks, Inc.
%   $Revision: 1.13.4.7 $  $Date: 2009/11/16 22:26:19 $

BI = [
    1       -183/64      37/12       -145/128
    0          0           0            0
    0       1500/371    -1000/159    1000/371
    0       -125/32       125/12     -375/64 
    0       9477/3392   -729/106    25515/6784
    0        -11/7        11/3        -55/28
    0         3/2         -4            5/2
    ];

s = (tinterp - t)/h;  
yinterp = y(:,ones(size(tinterp))) + f*(h*BI)*cumprod([s;s;s;s]);

end  


