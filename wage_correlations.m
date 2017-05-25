function [c,ceq]=wage_correlations(x)
%{  
    This function calculates the correlation of wage shocks between 
    partners in the household. It delivers:

    c       a vector of nonlinear inequalities c(x)<=0
    ceq     a vector of nonlinear equalities ceq(x)=0

    Alexandros Theloudis, LISER & UCL
    Email: a.theloudis@gmail.com

    -----------------------------------------------------------------------
%}

%   Declare parameters:
vH   = x(1) ;  % var men's permanent shock (Husband) 
uH   = x(2) ;  % var men's transitory shock
vW   = x(3) ;  % var women's permanent shock (Wife)
uW   = x(4) ;  % var women's transitory shock
vHW  = x(5) ;  % covar men's-women's permanent shocks
uHW  = x(6) ;  % covar men's-women's transitory shocks
rvHW = x(7) ;  % correlation: permanent shocks
ruHW = x(8) ;  % correlation: transitory shocks


%%  1.  NONLINEAR INEQUALITY CONSTRAINTS
%   Declare nonlinear inequality constraints c(x)<=0
%   -----------------------------------------------------------------------

c = [] ;


%%  2.  NONLINEAR EQUALITY CONSTRAINTS
%   Declare nonlinear equality constraints ceq(x)=0
%   -----------------------------------------------------------------------

%   Declare vector of nonlinear equality constraint:
%   (1 correlation of perm shocks and 1 correlation of trans shocks).
%   (A better code would have the model pass the number of constraints 
%   automatically to the function here, instead of hard-coding).
ceq = zeros(2,1) ;

%   Populate the vector of nonlinear constraints:
%   Correlation of permanent shocks
ceq(1) = rvHW - vHW / sqrt(vH*vW) ;
%   Correlation of transitory shocks
ceq(2) = ruHW - uHW / sqrt(uH*uW);

end