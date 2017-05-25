
/*******************************************************************************

Filename: 	LISER_MOMENTS_Stata_Shell.do
Author: 	Alexandros Theloudis (a.theloudis@gmail.com); LISER and UCL
Course:     Moments-Based Econometrics; 
			Luxembourg Institute for Socio-Economic Research
Date: 		March 2017

This file: 
(1)		Imports pre-processed wage and earnings data from the 
		US Panel Study of Income Dynamucs in years 1999-2011
(2)		Carries out multiple interactive GMM estimations of joint 
		wage dynamics of men and women (permanent-transitory process)

There is no copyright accompanying this file. Feel free to replicate, post 
online, or otherwise use as you wish. Please credit the author when you do so.

*******************************************************************************/

*	Initial statements:
clear
clear matrix
set more off
cap log close
cap version 14.0

*	Select computer for running code:
*	'=1' for Windows; '=2' for Unix.
global computer = 1

*	Set the path STATA will run from:
if $computer == 1 global DATAdir = "C:\"		/****** user provide input here! ******/
if $computer == 2 global DATAdir = "/" 			/****** user provide input here! ******/

*	Read PSID data: 
*	These are panel data (covering years 1999-2011 biennially) on wages and 
*	earnings of married men and women in the household. A household is the 
*	unit of observation over multiple years.
*
*	VARIABLES observed:
*	hh_id:	household id (cluster)
*	year:	year of observation (retrospective)
*	ageH: 	age of male household member
*	ageW: 	age of female household member
*	drwH:   residual growth in male wages (*see below for details*)
*	drwW:   residual growth in female wages
*	dreH:   residual growth in male earnings
*	dreW:   residual growth in female earnings
*	i_drXX: binary indicator ==1 if corresponding drXX above is non-missing
*			where XX = {wH,wW,eH,eW}
qui cd "$DATAdir"
if $computer == 1 import excel "$DATAdir\data\PSID_97_11_selected.xlsx", sheet("Sheet1") firstrow
if $computer == 2 import excel "$DATAdir/data/PSID_97_11_selected.xlsx", sheet("Sheet1") firstrow

*	Declare longitudinal nature of data:
*	hh_id is the household id (cluster); year (biennially) is the time variable.
tsset hh_id year, delta(2)

*	MODEL for WAGE DYNAMICS:
*		lnWageH(t) 					= X(t)'b + lnPermWageH(t) + uH(t)
*		lnPermWageH(t) 				= lnPermWageH(t-1) + vH(t)
*		lnWageH(t) - lnWageH(t-1) 	= (X(t)'-X(t-1)')b + vH(t) + uH(t) - uH(t-1)
*	where
*	-- lnWageH(t) is log real wage of the male household member at time t; 
*   -- X(t)'b captures observable effects including age; 
* 	-- lnPermWageH(t) is an AR(1) permanent wage component; 
*	-- vH(t) is the permanent shock at time t; 
*	-- uH(t) is a mean-reverting transitory shock. 
* 	I call 
*		drwH = lnWageH(t) - lnWageH(t-1) - (X(t)'-X(t-1)')b 
*	'residual growth in male wages' and similarly for the other variables.
*	Standard recent references include:
*	Meghir and Pistaferri (2004); Blundell et al. (2008). 

*	Generate lags:
*	Wage and earnings data are already first-differenced (i.e. expressed as 
*	growth over time or change in log; see 'MODEL' above for definitions).
*	Identification of the wage process requires fitting the covariance matrix 
*	of wage growth. Lags are needed in order to estimate the covariance between
*	wage growth in consecutive periods. The dataset appears artificially 
*	balanced. Whenever i_drwH == 0 (by construction == i_drwW == i_dreH == i_dreW), 
*	the corresponding observation is missing from the data.
cap drop if i_drwH == 0
foreach var of varlist drwH drwW {	
	qui gen L`var' = L.`var'
}
*

/*******************************************************************************
INTERACTIVE GMM ESTIMATION OF WAGE DYNAMICS
Stationarity of the wage process is imposed throughout
*******************************************************************************/

*	Estimation 1 - Equally Weighted GMM
*	Note the sample size and problems in convergence.
#delimit;
gmm 	(drwH^2 - {vH} - 2*{uH})		/* Var(drwH)(t) = vH(t) + uH(t) + uH(t-1) 			*/
		(drwH*LdrwH + {uH}) 			/* Cov(drwH(t),drwH(t-1)) = -uH(t-1)      			*/
		(drwW^2 - {vW} - 2*{uW})   		/* Var(drwW)(t) = vW(t) + uW(t) + uW(t-1) 			*/
		(drwW*LdrwW + {uW}) 			/* Cov(drwW(t),drwW(t-1)) = -uW(t-1)      			*/
		(drwH*drwW - {vHW} - 2*{uHW})   /* Cov(drwH,drwW)(t) = vHW(t) + uHW(t) + uHW(t-1) 	*/
		(LdrwH*drwW + {uHW})            /* Cov(drwH(t-1),drwW(t)) = -uHW(t-1)) 				*/
		(drwH*LdrwW + {uHW}), 			/* Cov(drwH(t),drwW(t-1)) = -uHW(t-1)) 				*/
		onestep winitial(identity) ;    /* identity weighting matrix */
#delimit cr

*	Estimation 2 - Equally Weighted GMM
*	Note the sample size. Am I using the right data now?
#delimit;
gmm 	(drwH^2 - {vH} - 2*{uH})
		(drwH*LdrwH + {uH})
		(drwW^2 - {vW} - 2*{uW})
		(drwW*LdrwW + {uW})
		(drwH*drwW - {vHW} - 2*{uHW})
		(LdrwH*drwW + {uHW})
		(drwH*LdrwW + {uHW}),
		onestep winitial(identity) nocommonesample ; /* does not restrict each moment to be evaluated on the same observations */
#delimit cr

*	Estimation 3 - Equally Weighted GMM without Instruments
*	Why STATA requires instruments in the first place?
*	Note STATA returns an error here.
#delimit;
capture 
gmm 	(drwH^2 - {vH} - 2*{uH})
		(drwH*LdrwH + {uH})
		(drwW^2 - {vW} - 2*{uW})
		(drwW*LdrwW + {uW})
		(drwH*drwW - {vHW} - 2*{uHW})
		(LdrwH*drwW + {uHW})
		(drwH*LdrwW + {uHW}),
		onestep winitial(identity) nocommonesample 
		instruments(,nocons) ;  /* for STATA, no instruments means no moments */
#delimit cr

*	Estimation 4 - Equally Weighted GMM; specify starting point
#delimit;
gmm 	(drwH^2 - {vH} - 2*{uH})
		(drwH*LdrwH + {uH})
		(drwW^2 - {vW} - 2*{uW})
		(drwW*LdrwW + {uW})
		(drwH*drwW - {vHW} - 2*{uHW})
		(LdrwH*drwW + {uHW})
		(drwH*LdrwW + {uHW}),
		onestep winitial(identity) nocommonesample 
		from(vH 0.05 uH 0.03 vW 0.05 uW 0.03 vHW 0.01 uHW 0.01) ;
#delimit cr

*	Estimation 5 - Equally Weighted GMM; request fast numerical derivatives
*	Where, what, and why are derivatives used in the first place?
#delimit;
capture
gmm 	(drwH^2 - {vH} - 2*{uH})
		(drwH*LdrwH + {uH})
		(drwW^2 - {vW} - 2*{uW})
		(drwW*LdrwW + {uW})
		(drwH*drwW - {vHW} - 2*{uHW})
		(LdrwH*drwW + {uHW})
		(drwH*LdrwW + {uHW}),
		onestep winitial(identity) nocommonesample quickd	/* quickd: calculates 'fast' derivatives; requires version 14. Note: these derivatives may fail in complicated nonlinear models */
		from(vH 0.05 uH 0.03 vW 0.05 uW 0.03 vHW 0.01 uHW 0.01) ;
#delimit cr

*	Estimation 6 - Equally Weighted GMM; supply exact analytical derivatives
*	Where, what, and why are derivatives used in the first place?
#delimit;
gmm 	(drwH^2 - {vH} - 2*{uH})
		(drwH*LdrwH + {uH})
		(drwW^2 - {vW} - 2*{uW})
		(drwW*LdrwW + {uW})
		(drwH*drwW - {vHW} - 2*{uHW}) 
		(LdrwH*drwW + {uHW})
		(drwH*LdrwW + {uHW}),
		onestep winitial(identity) nocommonesample 
		from(vH 0.05 uH 0.03 vW 0.05 uW 0.03 vHW 0.01 uHW 0.01)
		/* the Jacobian matrix of the model follows: */
		derivative(1/vH = -1)
		derivative(1/uH = -2) derivative(2/uH = 1)
		derivative(3/vW = -1)
		derivative(3/uW = -2) derivative(4/uW = 1) 
		derivative(5/vHW = -1)
		derivative(5/uHW = -2) derivative(6/uHW = 1) derivative(7/uHW = 1)  ;
#delimit cr

*	Estimation 7 - Optimal GMM without analytical derivatives
*	Note that Optimal GMM requires a two-stage estimation. 
#delimit;
gmm 	(drwH^2 - {vH} - 2*{uH})
		(drwH*LdrwH + {uH})
		(drwW^2 - {vW} - 2*{uW})
		(drwW*LdrwW + {uW})
		(drwH*drwW - {vHW} - 2*{uHW}) 
		(LdrwH*drwW + {uHW})
		(drwH*LdrwW + {uHW}),
		/* 'twostep' requests a two-stage estimator */
		twostep 
		/* winitial(identity) specifies the first-stage weighting matrix */
		winitial(identity) 
		/* vce(unadjusted) specifies type of covar matrix of parameters; here it is Optimal GMM */
		vce(unadjusted) nocommonesample ;
#delimit cr

*	Estimation 8 - Two-stage GMM with 'robust' standard errors
*	Note that 'robust' here accounts for errors that are independent but not 
*	necessarily identically distributed.
#delimit;
gmm 	(drwH^2 - {vH} - 2*{uH})
		(drwH*LdrwH + {uH})
		(drwW^2 - {vW} - 2*{uW})
		(drwW*LdrwW + {uW})
		(drwH*drwW - {vHW} - 2*{uHW}) 
		(LdrwH*drwW + {uHW})
		(drwH*LdrwW + {uHW}),
		twostep winitial(identity) wmatrix(robust) nocommonesample  ;
#delimit cr

*	Estimation 9 - Two-stage GMM with robust clustered standard errors
*	Note that 'cluster' accounts for arbitrary correlation among observations 
*	within clusters identified by a cluster variable (here hh_id).
#delimit;
gmm 	(drwH^2 - {vH} - 2*{uH})
		(drwH*LdrwH + {uH})
		(drwW^2 - {vW} - 2*{uW})
		(drwW*LdrwW + {uW})
		(drwH*drwW - {vHW} - 2*{uHW}) 
		(LdrwH*drwW + {uHW})
		(drwH*LdrwW + {uHW}),
		twostep winitial(identity) vce(cluster hh_id) nocommonesample  ;
#delimit cr

*	Estimation 10 - Two-stage GMM with robust clustered standard errors and 
*	exact analytical derivatives
#delimit;
gmm 	(drwH^2 - {vH} - 2*{uH})
		(drwH*LdrwH + {uH})
		(drwW^2 - {vW} - 2*{uW})
		(drwW*LdrwW + {uW})
		(drwH*drwW - {vHW} - 2*{uHW}) 
		(LdrwH*drwW + {uHW})
		(drwH*LdrwW + {uHW}),
		twostep winitial(identity) vce(cluster hh_id) nocommonesample 
		from(vH 0.05 uH 0.03 vW 0.05 uW 0.03 vHW 0.01 uHW 0.01)
		derivative(1/vH = -1)
		derivative(1/uH = -2) derivative(2/uH = 1)
		derivative(3/vW = -1)
		derivative(3/uW = -2) derivative(4/uW = 1) 
		derivative(5/vHW = -1)
		derivative(5/uHW = -2) derivative(6/uHW = 1) derivative(7/uHW = 1) ;
#delimit cr

*	Estimation 11 -  Equally Weighted GMM with bootstrap standard erros
*	Note that exact analytical derivatives are still supplied even though the 
*	standard errors no longer require these derivatives. Why?
preserve
#delimit;
gmm 	(drwH^2 - {vH} - 2*{uH})
		(drwH*LdrwH + {uH})
		(drwW^2 - {vW} - 2*{uW})
		(drwW*LdrwW + {uW})
		(drwH*drwW - {vHW} - 2*{uHW}) 
		(LdrwH*drwW + {uHW})
		(drwH*LdrwW + {uHW}),
		onestep winitial(identity) vce(bootstrap, reps(20) seed(071184)) nocommonesample 
		from(vH 0.05 uH 0.03 vW 0.05 uW 0.03 vHW 0.01 uHW 0.01)		
		derivative(1/vH = -1)
		derivative(1/uH = -2) derivative(2/uH = 1)
		derivative(3/vW = -1)
		derivative(3/uW = -2) derivative(4/uW = 1) 
		derivative(5/vHW = -1)
		derivative(5/uHW = -2) derivative(6/uHW = 1) derivative(7/uHW = 1) ;
#delimit cr
restore

*	Estimation 12 -  Equally Weighted GMM; block bootstrap standard errors 
*	(correct standard errors due to the longitudinal nature of data;
*	caveat: we are using pre-estimated residuals without accounting for this).
*	Reference for block bootstrap: Horowitz (2001).
preserve
#delimit;
gmm 	(drwH^2 - {vH} - 2*{uH})
		(drwH*LdrwH + {uH})
		(drwW^2 - {vW} - 2*{uW})
		(drwW*LdrwW + {uW})
		(drwH*drwW - {vHW} - 2*{uHW}) 
		(LdrwH*drwW + {uHW})
		(drwH*LdrwW + {uHW}),
		onestep winitial(identity) vce(bootstrap, reps(20) seed(071184) cluster(hh_id) idcluster(hh_id_bs) group(hh_id)) nocommonesample 
		from(vH 0.05 uH 0.03 vW 0.05 uW 0.03 vHW 0.01 uHW 0.01)		
		derivative(1/vH = -1)
		derivative(1/uH = -2) derivative(2/uH = 1)
		derivative(3/vW = -1)
		derivative(3/uW = -2) derivative(4/uW = 1) 
		derivative(5/vHW = -1)
		derivative(5/uHW = -2) derivative(6/uHW = 1) derivative(7/uHW = 1) ;
#delimit cr
restore

***** End of do file *****
