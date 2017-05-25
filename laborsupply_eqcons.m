function [Aeq,beq] = laborsupply_eqcons(x,mratios)
%{  
    This function constructs the equality constraint implied by 
    symmetry of the Frisch substitution matrix. This is an important 
    technical point that interested subjects can read on in Phlips (1974) 
    or other advanced textbooks in microeconomics. 
   
    The function delivers matrix [Aeq] and vector [beq] which are used as 
    linear equality constraints in the fmincon procedure.

    The equality constraint is of the form 
    Aeq * x = beq
    where Aeq is an MxN matrix (M linear equality constraints, N
    parameters) and beq is an Mx1 vector.

    Alexandros Theloudis, LISER & UCL
    Email: a.theloudis@gmail.com

    -----------------------------------------------------------------------
%}


%%  1.  SYMMETRY OF FRISCH SUBSTITUTION MATRIX
%   Define constraints (one constraint only in this case).
%   -----------------------------------------------------------------------

%   Frisch symmetry in labor supply cross-elasticities: 
%   eta_h2_w1 - eta_h1_w2 * E[W1*H1/W2*H2] = 0
Aeq1    = zeros(1,length(x)) ;
Aeq1(3) = 1 ;
Aeq1(2) = -mratios(1) ;


%%  2.  DELIVER
%   Stack constraints together (one constraint only in this case).
%   -----------------------------------------------------------------------

Aeq = Aeq1 ;
beq = zeros(size(Aeq,1),1) ;

end