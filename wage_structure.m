function f = wage_structure(vParams,mD,mI)
%{
    This function generates the moment conditions for the estimation of 
    the parameters of the wage process. It takes as arguments the vector of
    parameters [vParams], the data [mD], and the data indicators [mI].
 
    Note that the function calculates the reduced-form parameters for every
    iteration of the optimzation routine; this is redundant and the code 
    is therefore not the most efficient.

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
%   Initialize parameter vector and vectors of matching moments.
%   -----------------------------------------------------------------------

%   Parameter vector (second moments of wage shocks):
vH  = vParams(1) ;   % var men's permanent shock (Husband) 
uH  = vParams(2) ;   % var men's transitory shock
vW  = vParams(3) ;   % var women's permanent shock (Wife)
uW  = vParams(4) ;   % var women's transitory shock
vHW = vParams(5) ;   % covar men's-women's permanent shocks
uHW = vParams(6) ;   % covar men's-women's transitory shocks

%   Vectors to hold wage second moments:
V_wageH         = zeros(T,1) ;    %   VARIANCE of wage growth Husband  
C_wageHa        = zeros(T,1) ;    %   COVARIANCE of consecutive wage growths Husband
C_wageH_wageW   = zeros(T,1) ;    %   COVARIANCE of contemporaneous wage growths Husband-Wife
C_wageH_wageWa  = zeros(T,1) ;    %   COVARIANCE of consecutive wage growths Husband(t-1)-Wife(t)
C_wageH_wageWb  = zeros(T,1) ;    %   COVARIANCE of consecutive wage growths Husband(t)-Wife(t-1)
V_wageW         = zeros(T,1) ;    %   VARIANCE of wage growth Wife
C_wageWa        = zeros(T,1) ;    %   COVARIANCE of consecutive wage growths Wife


%%  2.  STRUCTURE OF MATCHING MOMENTS: INCLUDES REDUCED FORM AND STRUCTURAL PARAMETERS
%   Difference between empirical moments (reduced form parameters pi_hat) 
%   and structural counterpart g(theta_hat) must be 0:
%       pi = g(theta)
%   The code constructs the sample analog
%       pi_hat - g(theta_hat) = 0.
%   -----------------------------------------------------------------------

%   VARIANCE of wage growth husband:
%   Var(t)[DwH(i)] for periods t=2,3,4,5,6

%   Start the clock at t=2:
t = 2 ;

%   Start a loop that calculates Var(t)[DwH(i)] per t. Stop when t=T-1:
while t <= T-1
    
    %   For a given t, generate square of DwH(i) per observation. Note that
    %   this includes observations for which DwH(i) is missing (==0):
    wi = mD(t,:).*mD(t,:) ;
    
    %   For a given t, generate indicators such that indicator==1 means
    %   that DwH(i) is not missing for observation i:
    ni = mI(t,:).*mI(t,:) ;
    
    %   Get the mean of (DwH(i))^2 across i (thus the variance) only using 
    %   observations for which indicator==1
    V_wageH(t) = mean(wi(ni==1)) ;
    
    %   Move on to next period:
    t = t + 1 ;
end 

%   The theoretical counterpart of the empirical variance Var(t)[DwH(i)] 
%   is the sum of the variance of the permanent shock (vH) and twice the 
%   variance of the transitory shock (uH) -- this only holds under 
%   stationarity of the wage process. Thus the difference between the 
%   empirical variance V_wageH and its theoretical counterpart must be 0. 
%   This is one (vector of) moment(s) used in the estimation: 
V_wageH = V_wageH - vH - 2*uH ;


%   -----------------------------------------------------------------------
%   COVARIANCE of consecutive wage growths husband:
%   Cov[DwH(i)(t), DwH(i)(t-1)] for periods t=2,3,4,5,6,7

t = 2 ;
while t <= T
    wi = mD(t,:).*mD(t-1,:) ;
    ni = mI(t,:).*mI(t-1,:) ;
    C_wageHa(t) = mean(wi(ni==1)) ;
    t = t + 1 ;
end
C_wageHa = C_wageHa  + uH ;


%   -----------------------------------------------------------------------
%   COVARIANCE of contemporaneous wage growths husband-wife:
%   Cov[DwH(i)(t), DwW(i)(t)] for periods t=2,3,4,5,6

t = 2 ;
while t <= T-1
    wi = mD(t,:).*mD(T+t,:) ;
    ni = mI(t,:).*mI(T+t,:) ;
	C_wageH_wageW(t) = mean(wi(ni==1)) ;
	t = t + 1 ;
end
C_wageH_wageW = C_wageH_wageW - vHW - 2*uHW ;
    

%   -----------------------------------------------------------------------
%   COVARIANCE of consecutive wage growths husband-wife:
%   Cov[DwH(i)(t-1), DwW(i)(t)] for periods t=2,3,4,5,6,7

t = 2 ;
while t <= T
    wi = mD(t-1,:).*mD(T+t,:) ;
    ni = mI(t-1,:).*mI(T+t,:) ;
    C_wageH_wageWa(t) = mean(wi(ni==1)) ;
	t = t + 1 ;
end
C_wageH_wageWa = C_wageH_wageWa + uHW ;


%   -----------------------------------------------------------------------
%   COVARIANCE of consecutive wage growths husband-wife:
%   Cov[DwH(i)(t), DwW(i)(t-1)] for periods t=2,3,4,5,6,7

t = 2 ;
while t <= T 
    wi = mD(t,:).*mD(T+t-1,:) ;
    ni = mI(t,:).*mI(T+t-1,:) ;
	C_wageH_wageWb(t) = mean(wi(ni==1)) ; 
	t = t + 1 ;
end
C_wageH_wageWb = C_wageH_wageWb + uHW ; 


%   -----------------------------------------------------------------------
%   VARIANCE of wage growth wife:
%   Var(t)[DwW(i)] for periods t=2,3,4,5,6

t = 2 ;
while t <= T-1
    wi = mD(T+t,:).*mD(T+t,:) ;
    ni = mI(T+t,:).*mI(T+t,:) ;
    V_wageW(t) = mean(wi(ni==1)) ;
    t = t + 1 ;
end 
V_wageW = V_wageW  - vW - 2*uW ;
    

%   -----------------------------------------------------------------------
%   COVARIANCE of consecutive wage growths wife:
%   Cov[DwW(i)(t), DwW(i)(t-1)] for periods t=2,3,4,5,6,7

t = 2 ;
while t <= T
    wi = mD(T+t,:).*mD(T+t-1,:) ;
    ni = mI(T+t,:).*mI(T+t-1,:) ;
    C_wageWa(t) = mean(wi(ni==1)) ;
    t = t + 1 ;
end
C_wageWa = C_wageWa + uW ;


%%  3.  STACK MOMENTS TOGETHER
%   Stack vectors of moments together, remove 0 entries, deliver objective 
%   function (quadratic loss).
%   -----------------------------------------------------------------------

%   Stack together:
vMoms = [   V_wageH(2:T-1)  ; ...      
            C_wageHa(2:T)  ; ...
            C_wageH_wageW(2:T-1)  ; ...
            C_wageH_wageWa(2:T)  ; ...
            C_wageH_wageWb(2:T)  ; ...
            V_wageW(2:T-1)  ; ...
            C_wageWa(2:T)] ;

%   Objective function with identity as the weighting matrix to be returned 
%   to the optimizer:
f = vMoms' * eye(size(vMoms,1)) * vMoms ;