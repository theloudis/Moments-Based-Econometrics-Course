function f = laborsupply_structure(vParams,wagemoms,mD,mI)
%{  
    This function generates the moment conditions for the estimation of 
    the parameters of labor supply. It takes as arguments the vector of
    parameters [vParams], the vector of pre-estimated wage moments 
    [wagemoms], the wage and earnings data [mD], and the data indicators [mI].

    The function is commented throughout.

    There is no copyright accompanying this code. Feel free to replicate, 
    post online, or otherwise use as you wish. Please credit the author 
    when you do so. 

    Alexandros Theloudis, LISER & UCL
    Email: a.theloudis@gmail.com

    -----------------------------------------------------------------------
%}

%   Initial statements:
global T ;


%%  1.  INITIALIZE OBJECTS
%   Initialize parameters and vectors to hold orthogonality conditions.
%   -----------------------------------------------------------------------

%   Parameter vector (structural labor supply parameters):
eta_h1_w1 = vParams(1) ; % men's own-wage labor supply elasticity
eta_h1_w2 = vParams(2) ; % men's cross-wage labor supply elasticity
eta_h2_w1 = vParams(3) ; % women's cross-wage labor supply elasticity
eta_h2_w2 = vParams(4) ; % women's own-wage labor supply elasticity

%   Declare pre-estimated first-stage wage moments (note that these moments 
%   have not been corrected for wage measurement error -- this contributes 
%   to estimating negative labor supply elasticities from the response of
%   earnings to transitory shocks alone):
uH  = wagemoms(2) ; % var men's transitory shock
uW  = wagemoms(4) ; % var women's transitory shock
uHW = wagemoms(6) ; % covar men's-women's transitory shocks

%   Vectors to hold (part of the) covariance matrix of earnings and wages:     
C_wageH_yHa = zeros(T,1) ; % COVARIANCE of consecutive wage (t-1) and earnings (t) growth Husband
C_wageH_yHb = zeros(T,1) ; % COVARIANCE of consecutive wage (t) and earnings (t-1) growth Husband
C_wageW_yHa = zeros(T,1) ; % COVARIANCE of female wage (t-1) and male earnings (t) growth
C_wageW_yHb = zeros(T,1) ; % COVARIANCE of female wage (t) and male earnings (t-1) growth
C_wageH_yWa = zeros(T,1) ; % COVARIANCE of male wage (t-1) and female earnings (t) growth
C_wageH_yWb = zeros(T,1) ; % COVARIANCE of male wage (t) and female earnings (t-1) growth    
C_wageW_yWa = zeros(T,1) ; % COVARIANCE of consecutive wage (t-1) and earnings (t) growth Wife
C_wageW_yWb = zeros(T,1) ; % COVARIANCE of consecutive wage (t) and earnings (t-1) growth Wife
C_yH_yHa    = zeros(T,1) ; % COVARIANCE of consecutive earnings (t-1,t) growth Husband
C_yH_yWa    = zeros(T,1) ; % COVARIANCE of male earnings (t-1) and female earnings (t) growth
C_yH_yWb    = zeros(T,1) ; % COVARIANCE of male earnings (t) and female earnings (t-1) growth  
C_yW_yWa    = zeros(T,1) ; % COVARIANCE of consecutive earnings (t-1,t) growth Wife


%%  2.  BUILD ORTHOGONALITY CONDITIONS
%   Construct vector of orthogonality conditions implied by the structural
%   model such that
%   
%   E[h(w_i,theta)] = 0
%   
%   This part of the code constructs vector h(.) as function of the data 
%   and structural parameter theta.
%   -----------------------------------------------------------------------

%   COVARIANCE of consecutive wage (t-1) and earnings (t) growth Husband:
%   Cov[DeH(t-1),DyH(t)] for periods t=2,3,4,5,6,7

%   Start the clock at t=2:
t = 2 ;

%   Start a loop that calculates Cov[DwH(t-1),DyH(t)] per t. Stop when t=T.
while t <= T 
    
    %   For a given t, request position in vector C_wageH_yHa whose 
    %   corresponding orthogonality condition we'll now build:
	C_wageH_yHa(t) = ( ...                              
        (mD(t-1,:).*mD(2*T+t,:)) ...                    % data
        + (1+eta_h1_w1) * uH + eta_h1_w2 * uHW) ...     % structure
        * (mI(t-1,:).*mI(2*T+t,:))' ;                   % adjustment for unbalanced panel / missing observation
    
    %   Sample analog of population moment requires that we adjust for
    %   (divide by) sample size:
    C_wageH_yHa(t) = ...
        C_wageH_yHa(t) / (mI(t-1,:)*mI(2*T+t,:)') ; 
    
    %   Move on to next period:
	t = t + 1 ;
end


%   -----------------------------------------------------------------------
%   COVARIANCE of consecutive wage (t) and earnings (t-1) growth Husband:
%   Cov[DwH(t),DyH(t-1)] for periods t=2,3,4,5,6,7

t = 2 ;
while t <= T 
    C_wageH_yHb(t) = ((mD(t,:).*mD(2*T+t-1,:)) ...
        + (1+eta_h1_w1) * uH + eta_h1_w2 * uHW) ...
        * (mI(t,:).*mI(2*T+t-1,:))' ;
    C_wageH_yHb(t) = C_wageH_yHb(t) / (mI(t,:)*mI(2*T+t-1,:)') ;
	t = t + 1 ;
end

   
%   -----------------------------------------------------------------------
%   COVARIANCE of female wage (t-1) and male earnings (t) growth
%   Cov[DwW(t-1),DyH(t)] for periods t=2,3,4,5,6,7

t = 2 ;
while t <= T 
	C_wageW_yHa(t) = ((mD(T+t-1,:).*mD(2*T+t,:)) ...
        + (1+eta_h1_w1) * uHW + eta_h1_w2 * uW) ...
        * (mI(T+t-1,:).*mI(2*T+t,:))' ;
    C_wageW_yHa(t) = C_wageW_yHa(t) / (mI(T+t-1,:)*mI(2*T+t,:)') ;
	t = t + 1 ;
end


%   -----------------------------------------------------------------------
%   COVARIANCE of female wage (t) and male earnings (t-1) growth
%   Cov[DwW(t),DyH(t-1)] for periods t=2,3,4,5,6,7

t = 2 ;
while t <= T
    C_wageW_yHb(t) = ((mD(T+t,:).*mD(2*T+t-1,:)) ...
        + (1+eta_h1_w1) * uHW + eta_h1_w2 * uW) ...
        * (mI(T+t,:).*mI(2*T+t-1,:))' ;
    C_wageW_yHb(t) = C_wageW_yHb(t) / (mI(T+t,:)*mI(2*T+t-1,:)') ;
	t = t + 1 ;
end


%   -----------------------------------------------------------------------
%   COVARIANCE of male wage (t-1) and female earnings (t) growth
%   Cov[DwH(t-1),DyW(t)] for periods t=2,3,4,5,6,7

t = 2 ;
while t <= T 
	C_wageH_yWa(t) = ((mD(t-1,:).*mD(3*T+t,:)) ...
        + eta_h2_w1 * uH + (1+eta_h2_w2) * uHW) ...
        * (mI(t-1,:).*mI(3*T+t,:))' ; 
    C_wageH_yWa(t) = C_wageH_yWa(t) / (mI(t-1,:)*mI(3*T+t,:)') ;
	t = t + 1 ;
end


%   -----------------------------------------------------------------------
%   COVARIANCE of male wage (t) and female earnings (t-1) growth
%   Cov[DwH(t),DyW(t-1)] for periods t=2,3,4,5,6,7

t = 2 ;
while t <= T 
	C_wageH_yWb(t) = ((mD(t,:).*mD(3*T+t-1,:)) ...
        + eta_h2_w1 * uH + (1+eta_h2_w2) * uHW) ...
        * (mI(t,:).*mI(3*T+t-1,:))' ;
    C_wageH_yWb(t) = C_wageH_yWb(t) / (mI(t,:)*mI(3*T+t-1,:)') ;
	t = t + 1 ;
end


%   -----------------------------------------------------------------------
%   COVARIANCE of consecutive wage (t-1) and earnings (t) growth Wife   
%   Cov[DwW(t-1),DyW(t)] for periods t=2,3,4,5,6,7

t = 2 ;
while t <= T 
	C_wageW_yWa(t) = ((mD(T+t-1,:).*mD(3*T+t,:)) ...
        + eta_h2_w1 * uHW + (1+eta_h2_w2) * uW) ...
        * (mI(T+t-1,:).*mI(3*T+t,:))' ;
    C_wageW_yWa(t) = C_wageW_yWa(t) / (mI(T+t-1,:)*mI(3*T+t,:)') ;
	t = t + 1 ;
end


%   -----------------------------------------------------------------------
%   COVARIANCE of consecutive wage (t) and earnings (t-1) growth Wife   
%   Cov[DwW(t),DyW(t-1)] for periods t=2,3,4,5,6,7

t = 2 ;
while t <= T 
	C_wageW_yWb(t) = ((mD(T+t,:).*mD(3*T+t-1,:)) ...
        + eta_h2_w1 * uHW + (1+eta_h2_w2) * uW) ...
        * (mI(T+t,:).*mI(3*T+t-1,:))' ;
    C_wageW_yWb(t) = C_wageW_yWb(t) / (mI(T+t,:)*mI(3*T+t-1,:)') ;
	t = t + 1 ;
end


%   -----------------------------------------------------------------------
%   COVARIANCE of consecutive earnings (t-1,t) growth Husband
%   Cov[DyH(t-1),DyH(t)] for periods t=2,3,4,5,6,7

t = 2 ;
while t <= T 
	C_yH_yHa(t) = ((mD(2*T+t-1,:).*mD(2*T+t,:)) ...
        + (1 + eta_h1_w1)^2 * uH ...
        + 2*eta_h1_w2*(1 + eta_h1_w1) * uHW ...
        + eta_h1_w2^2 * uW) ...
        * (mI(2*T+t-1,:).*mI(2*T+t,:))' ;
    C_yH_yHa(t) = C_yH_yHa(t) / (mI(2*T+t-1,:)*mI(2*T+t,:)') ;
	t = t + 1 ;
end

 
%   -----------------------------------------------------------------------
%   COVARIANCE of male earnings (t-1) and female earnings (t) growth
%   Cov[DyH(t-1),DyW(t)] for periods t=2,3,4,5,6,7

t = 2 ;
while t <= T
	C_yH_yWa(t) = ((mD(2*T+t-1,:).*mD(3*T+t,:)) ...
        + eta_h2_w1 * (1 + eta_h1_w1) * uH ...
        + (1 + eta_h2_w2)*(1 + eta_h1_w1) * uHW ...
        + eta_h1_w2*eta_h2_w1 * uHW ...
        + eta_h1_w2 * (1 + eta_h2_w2) * uW) ...    
        * (mI(2*T+t-1,:).*mI(3*T+t,:))' ; 
    C_yH_yWa(t) = C_yH_yWa(t) / (mI(2*T+t-1,:)*mI(3*T+t,:)') ;
	t = t + 1 ;
end


%   -----------------------------------------------------------------------
%   COVARIANCE of male earnings (t) and female earnings (t-1) growth
%   Cov[DyH(t),DyW(t-1)] for periods t=2,3,4,5,6,7

t = 2 ;
while t <= T
	C_yH_yWb(t) = ((mD(2*T+t,:).*mD(3*T+t-1,:)) ...
        + eta_h2_w1 * (1 + eta_h1_w1) * uH ...
        + (1 + eta_h1_w1)*(1 + eta_h2_w2) * uHW ...
        + eta_h1_w2*eta_h2_w1 * uHW ...
        + eta_h1_w2 * (1 + eta_h2_w2) * uW) ... 
        * (mI(2*T+t,:).*mI(3*T+t-1,:))' ;
    C_yH_yWb(t) = C_yH_yWb(t) / (mI(2*T+t,:)*mI(3*T+t-1,:)') ;
	t = t + 1 ;
end


%   -----------------------------------------------------------------------
%   COVARIANCE of consecutive earnings (t-1,t) growth Wife
%   Cov[DyW(t-1),DyW(t)] for periods t=2,3,4,5,6,7

t = 2 ;
while t <= T 
	C_yW_yWa(t) = ((mD(3*T+t-1,:).*mD(3*T+t,:)) ...
        + eta_h2_w1^2 * uH ...
        + 2*eta_h2_w1 * (1 + eta_h2_w2) * uHW ...
        + (1 + eta_h2_w2)^2 * uW) ...
        * (mI(3*T+t-1,:).*mI(3*T+t,:))' ;
    C_yW_yWa(t) = C_yW_yWa(t) / (mI(3*T+t-1,:)*mI(3*T+t,:)') ;
	t = t + 1 ;
end


%%  3.  DELIVER
%   Stack vectors of moments together, remove 0 entries*, deliver objective 
%   function (quadratic loss) to pass to optimizer.
%   *not needed when using identity matrix
%   -----------------------------------------------------------------------
vMoms = [   C_wageH_yHa ; ...
            C_wageH_yHb ; ...
            C_wageW_yHa ; ...
            C_wageW_yHb ; ...
            C_wageH_yWa ; ...
            C_wageH_yWb ; ...
            C_wageW_yWa ; ...
            C_wageW_yWb ; ...        
            C_yH_yHa    ; ...
            C_yH_yWa    ; ...
            C_yH_yWb    ; ...
            C_yW_yWa  ] ;
        
%   Objective function criterion:
f = vMoms' * eye(size(vMoms,1)) * vMoms ;