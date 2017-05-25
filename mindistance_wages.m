function [what,wfval,wflag] = mindistance_wages(mD,mInd)
%{  
    This routine implements the Minimum Distance estimation of the 
    parameters of the wage process. It accepts data matrix [mD] and data
    indicators [mInd] and returns a vector of MD point estimates [what],
    the value of the objective function at the optimum (wfval) and the
    optimization exit flag (wflag). Caution if the exit flag is not 1.
 
	MODEL for WAGE DYNAMICS:
		lnWageH(t)                          = X(t)'b + lnPermWageH(t) + uH(t)
		lnPermWageH(t)                      = lnPermWageH(t-1) + vH(t)
		drwH = lnWageH(t) - lnWageH(t-1) 	= (X(t)'-X(t-1)')b + vH(t) + uH(t) - uH(t-1)
	where
	-- lnWageH(t) is log real wage of the male household member at time t; 
    -- X(t)'b captures observable effects including age; 
 	-- lnPermWageH(t) is an AR(1) permanent wage component; 
	-- vH(t) is the permanent shock at time t; 
	-- uH(t) is a mean-reverting transitory shock at time t. 
 	Similarly for female wages. Standard recent references include:
	Meghir and Pistaferri (2004); Blundell et al. (2008).
   
    There is no copyright accompanying this code. Feel free to replicate, 
    post online, or otherwise use as you wish. Please credit the author 
    when you do so. 

    Alexandros Theloudis, LISER & UCL
    Email: a.theloudis@gmail.com

    -----------------------------------------------------------------------
%}

%%  1.  ESTIMATION ATTRIBUTES
%   Define estimation features and initial parameter vector.
%   -----------------------------------------------------------------------

%   Estimation features:
%   This part is technical but important. For our specific solver, fmincon, read more: https://uk.mathworks.com/help/optim/ug/fmincon.html#input_argument_d0e49295 
wage_opt_alg            = 'interior-point' ;    % which optimization algorithm; read more: https://uk.mathworks.com/help/optim/ug/choosing-the-algorithm.html
wage_opt_contol         = 1e-6 ;                % constraint tolerance
wage_opt_displ          = 'iter' ;              % display of progress
wage_opt_funtol         = 1e-6 ;                % function tolerance
wage_opt_maxfunevl      = 1e+5 ;                % max number of function evaluations; allow big number if no analytical derivatives and #var is big
wage_opt_maxiter        = 1000 ;                % max number of iterations

%   Declare initial parameter values (a total of 8 parameters):
x0 = [ 0.05     ; ... % vH    (1) var men's permanent shock (Husband)   
       0.03     ; ... % uH    (2) var men's transitory shock
       0.05     ; ... % vW    (3) var women's permanent shock (Wife)
       0.03     ; ... % uW    (4) var women's transitory shock
       0.01     ; ... % vHW   (5) covar men's-women's permanent shocks
       0.01     ; ... % uHW   (6) covar men's-women's transitory shocks
       0.1      ; ... % rvHW  (7) Pearson correlation men's-women's permanent shocks
       0.1 ]    ; ... % ruHW  (8) Pearson correlation men's-women's transitory shocks
  
%   Declare parameter bounds:
xb = [ 0.0 1.0  ; ... % vH    (1)
       0.0 1.0  ; ... % uH    (2)
       0.0 1.0  ; ... % vW    (3)
       0.0 1.0  ; ... % uW    (4)
      -0.1 0.1  ; ... % vHW   (5)
      -0.1 0.1  ; ... % uHW   (6)
      -1.0 1.0  ; ... % rvHW  (7)
      -1.0 1.0 ]; ... % ruHW  (8)

              
%%  2.  IMPLEMENT MINIMUM DISTANCE
%   Impose estimation features and implement MD estimation.
%   -----------------------------------------------------------------------

%   Define optimisation options set:
fminconopt = optimoptions('fmincon', ...
                          'Algorithm',wage_opt_alg, ...
                          'ConstraintTolerance',wage_opt_contol, ...
                          'Display',wage_opt_displ, ...
                          'FunctionTolerance',wage_opt_funtol, ...
                          'MaxFunctionEvaluations',wage_opt_maxfunevl, ...
                          'MaxIterations',wage_opt_maxiter) ;

%   Implement MD optimisation calling MATLAB's fmincon:
%   fmincon is a nonlinear programming solver, it searches for the minimum 
%   of a constrained nonlinear multivariable function f: R^q -> R^1
%   Read more: https://uk.mathworks.com/help/optim/ug/fmincon.html
[what,wfval,wflag] = fmincon(@(x) wage_structure(x,mD,mInd), ...    % here I pass the objective function to be minimized
                                x0, ...                             % here I pass the parameter starting point
                                [],[],[],[], ...                    % here I pass linear equality and inequality constraints on the parameters
                                xb(:,1),xb(:,2), ...                % here I pass parameter bounds
                             @(x) wage_correlations(x),...          % here I pass nonlinear constraints
                                fminconopt) ;                       % here I pass the optimzation options
end