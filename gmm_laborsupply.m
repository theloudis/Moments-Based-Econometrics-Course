function [lhat,lfval,lflag] = gmm_laborsupply(wagemoms,mD,mInd)
%{ 
    This routine implements the Generalized Method of Moments estimation 
    of the parameters of family labor supply using earnings and wage 
    moments. Identification of the structure of labor supply relies on
    earnings (a function of labor supply) responding to exogenous variation 
    in wages (wage shocks). The rourine accepts wage moments [wagemoms], 
    data matrix [mD] and data indicators [mInd] and returns a vector of 
    GMM point estimates [lhat], the value of the objective function at the 
    optimum (lfval) and the optimization exit flag (lflag). Caution if 
    the exit flag is not 1.

    MODEL OF EARNINGS DYNAMICS
        D(lnEarningsH(t)) = eta_hH_wH * D(uH(t)) + eta_hH_wW * D(uW(t)) + gEarningsH(vH,vW)
        D(lnEarningsH(t)) = lnEarningsH(t) - lnEarningsH(t-1) - (X(t)'-X(t-1)')b
		D(uH(t))          = uH(t) - uH(t-1)
        D(uW(t))          = uW(t) - uW(t-1)
	where
	-- lnEarningsH(t) is log real earnings of the male household member at time t; 
    -- uH(t) (uW(t)) is a mean-reverting transitory shock at time t;
    -- eta_hH_wH is the male labor supply elasticity (Frisch) with respect to his own wage
    -- eta_hH_wW is the male labor supply elasticity (Frisch) with respect to his wife's wage
    -- gEarningsH(vH,vW) is a function of male/female permanent shocks
    -- X(t)'b captures observable effects including age; 
 	Similarly for female earnings. Recent reference: Blundell et al. (2016).
   
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
lsupply_opt_alg       = 'interior-point' ;  % which optimization algorithm; read more: https://uk.mathworks.com/help/optim/ug/choosing-the-algorithm.html
lsupply_opt_contol    = 1e-6 ;              % constraint tolerance
lsupply_opt_displ     = 'iter' ;            % display of progress
lsupply_opt_funtol    = 1e-6 ;              % function tolerance
lsupply_opt_maxfunevl = 2e+4 ;              % max number of function evaluations; allow big number if no analytical derivatives and #var is big
lsupply_opt_maxiter   = 1000 ;              % max number of iterations

%   Declare initial parameter values (a total of 4 parameters):
x0 = [ 0.5    % eta_h1_w1 (1) men's own-wage labor supply elasticity
       0.1    % eta_h1_w2 (2) men's cross-wage labor supply elasticity
       0.1    % eta_h2_w1 (3) women's cross-wage labor supply elasticity
       0.8 ]; % eta_h2_w2 (4) women's own-wage labor supply elasticity   

%   Declare parameter bounds:
xb = [ -2.0 2.0 ; ...   % eta_h1_w1 (1)
       -2.0 2.0 ; ...   % eta_h1_w2 (2)
       -2.0 2.0 ; ...   % eta_h2_w1 (3)
       -2.0 2.0];       % eta_h2_w2 (4)


%%  2.  IMPLEMENT GMM ESTIMATION OF LABOR SUPPLY MODEL
%   Impose estimation features and implement GMM estimation.
%   -----------------------------------------------------------------------

%   Define linear equality constraints on 'added-worker' parameters:
%   Here I impose Frisch symmetry in labor supply cross-elasticities, 
%   meaning 
%       eta_h2_w1 - eta_h1_w2 * E[W1*H1/W2*H2] = 0.
%   This part is technical and I refer to Phlips (1974).
%   Hard-coded value E[W1*H1/W2*H2] = 2.0932 is taken from the PSID data.
%   (A better code would have the data pass this number automatically to
%   the function here, instead of hard-coding).
[Aeq,beq] = laborsupply_eqcons(x0,2.0932) ;

%   Define optimisation options set:
fminconopt = optimoptions('fmincon', ...
                          'Algorithm',lsupply_opt_alg, ...
                          'ConstraintTolerance',lsupply_opt_contol, ...
                          'Display',lsupply_opt_displ, ...
                          'FunctionTolerance',lsupply_opt_funtol, ...
                          'MaxFunctionEvaluations',lsupply_opt_maxfunevl, ...
                          'MaxIterations',lsupply_opt_maxiter) ;

%   Implement GMM optimisation calling MATLAB's fmincon:
%   fmincon is a nonlinear programming solver, it searches for the minimum 
%   of a constrained nonlinear multivariable function f: R^q -> R^1
%   Read more: https://uk.mathworks.com/help/optim/ug/fmincon.html
[lhat,lfval,lflag] = fmincon(@(x) laborsupply_structure(x,wagemoms,mD,mInd), ...    % here I pass the objective function to be minimized
                                x0,...                                              % here I pass the parameter starting point
                                [],[],Aeq,beq, ...                                  % here I pass linear equality and inequality constraints on the parameters
                                xb(:,1),xb(:,2),...                                 % here I pass parameter bounds
                                [],...                                              % here I pass nonlinear constraints
                                fminconopt) ;                                       % here I pass the optimzation options

end