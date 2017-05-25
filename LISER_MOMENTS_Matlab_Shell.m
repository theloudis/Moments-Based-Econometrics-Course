
%%  'Moments-Based Methods for Structural Estimation: Theory and Applications'
%   
%   Filename: 	LISER_MOMENTS_Matlab_Shell.do
%   Author: 	Alexandros Theloudis (a.theloudis@gmail.com); LISER and UCL
%   Course:     Moments-Based Econometrics; 
%               Luxembourg Institute for Socio-Economic Research
%   Date: 		March 2017

%   This file: 
%   (1)	Imports pre-processed wage and earnings data from the 
%       US Panel Study of Income Dynamucs in years 1999-2011.
%   (2)	Carries out a Minimum Distance estimation of the joint 
%       wage dynamics of men and women (permanent-transitory process)
%   (3) Carries out a GMM estimation of labor supply preferences using 
%       analytical equations for men's and women's labor supply implied 
%       by an intertemporal model for consumption and labor supply.
%       This part uses results for wages from (2).

%   This script calls basic MATLAB functions as well as functions available
%   in the optimization toolbox. Other functions that do not fall in the
%   above two categories are provided with this distribution.

%   There is no copyright accompanying this code. Feel free to replicate, 
%   post online, or otherwise use as you wish. Please credit the author 
%   when you do so. File more easily implemented if read alongside 
%   accompanying course notes.

%   Alexandros Theloudis, LISER & UCL
%   Email: a.theloudis@gmail.com
%
%   =======================================================================


%%  1.  INITIAL STATEMENTS
%   Set the program switches & attributes; read the data.
%   -----------------------------------------------------------------------

clear
clc

global  computer do_labor_supply do_slices ...
        T winDir macDir ;

%   Switches:
computer = 2 ;          % computer = { 1 PC; 2 Mac; }
do_labor_supply = 1 ;   % do_labor_supply = { 0 No; 1 Yes }
do_slices = 1 ;         % do_slices = { 0 No; 1 Yes }

%   Directories:
disp('User: have you provided your working directory?')
pause ;
winDir = 'C:\' ;        % ******* user provide input! *****
macDir = '/' ;          % ******* user provide input! *****
if computer == 1        % windows          
    cd(winDir) ;
elseif computer == 2    % mac 
    cd(macDir) ;
else                    % error
    error('Invalid selection of computer.')
end
addpath('data') ;

%   Read PSID data:
T = 7 ;                 % number of periods in data
[mData,mIndicators] = readData('PSID_97_11_selected.xlsx') ;


%%  2.  MINIMUM DISTANCE ESTIMATION OF FAMILY WAGE DYNAMICS
%   This part implements the MD estimation of the paramaters of the wage 
%   process. Features and switches for this estimation are set 
%   inside function 'mindistance_wages'.
%   -----------------------------------------------------------------------

%   Call MD estimation routine:
[vWageHat,wageFval,wageFlag] = mindistance_wages(mData,mIndicators) ;

%   Rearrange estimates for ease of illustration:
mWageHat = [vWageHat(1)/2 vWageHat(3)/2 vWageHat(5)/2 vWageHat(7); ...
            vWageHat(2)   vWageHat(4)   vWageHat(6)  vWageHat(8)];

        
%%  3.  GMM ESTIMATION OF LABOR SUPPLY PARAMETERS
%   This part implements the GMM estimation of the paramaters of the 
%   structural model. Features and switches for this estimation are set 
%   inside script 'gmm_laborsupply'.
%   -----------------------------------------------------------------------

%   Call GMM routine:
if do_labor_supply == 1
    [vModelHat,modelFval,modelFlag] = gmm_laborsupply(vWageHat,mData,mIndicators) ;
end


%%  4.  CHECK PERFORMANCE OF IDENTIFICATION AROUND OPTIMAL AND GLOBAL OPTIMUM
%   This part plots slices of the objective f around the optimal point.
%   -----------------------------------------------------------------------

%   Assesses behavior of function, estimation, identification:
if do_slices == 1
    [xvls,fvls,feas] = slices_model(vModelHat,vWageHat,mData,mIndicators,0.3) ;
end