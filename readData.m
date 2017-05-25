function [rmdata,rmind] = readData(dataname)
%{
    This function reads the PSID data (dataname) and returns 2 matrices: 
    [rmdata], which contains the artificially balanced panel data, and 
    [rmind], which contains indicators for missing data. 

    There is no copyright accompanying this code. Feel free to replicate, 
    post online, or otherwise use as you wish. Please credit the author 
    when you do so. 

    Alexandros Theloudis, LISER & UCL
    Email: a.theloudis@gmail.com

    -----------------------------------------------------------------------
%}

%   Initial statements:
global computer T winDir macDir ;
if      computer == 1 % Windows          
    cd(winDir) ;
elseif  computer == 2 % Mac 
    cd(macDir) ;
end

%   Read data:
mdata = xlsread(dataname) ;

%	VARIABLES included:
%	hh_id:	household id (cluster)
%	year:	year of observation (retrospective)
%	ageH: 	age of male household member
%	ageW: 	age of female household member
%	drwH:   residual growth in male wages 
%	drwW:   residual growth in female wages
%	dreH:   residual growth in male earnings
%	dreW:   residual growth in female earnings
%	i_drXX: binary indicator ==1 if corresponding drXX above is non-missing
%			where XX = {wH,wW,eH,eW}

%   Declare variable position in dataset:
cl_drwH    = 5  ;
cl_i_drwH  = 6  ;
cl_drwW    = 7  ;
cl_i_drwW  = 8  ;
cl_dreH    = 9  ;
cl_i_dreH  = 10 ;
cl_dreW    = 11 ;
cl_i_dreW  = 12 ;

%   Obtain sample sizes and reshape data into wide format:
%   Now each column of length 28=4*T is one observation (household; hh_id) 
%   with T=7 waves of male wage data, followed by T waves of female wage 
%   data, followed by T waves of male earnings data, followed by T waves
%   of female earnings data.
l = size(mdata) ;
n = l(1)/T ;
rmdata = vertcat(   reshape(mdata(:,cl_drwH),T,n),      ...
                    reshape(mdata(:,cl_drwW),T,n),      ...  
                    reshape(mdata(:,cl_dreH),T,n),      ...   
                    reshape(mdata(:,cl_dreW),T,n));
rmind  = vertcat(   reshape(mdata(:,cl_i_drwH),T,n),    ... 
                    reshape(mdata(:,cl_i_drwW),T,n),    ...
                    reshape(mdata(:,cl_i_dreH),T,n),    ...
                    reshape(mdata(:,cl_i_dreW),T,n));
end