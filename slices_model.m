function [xvalues,fvalues,feasibl] = slices_model(phat,what,mD,mInd,fineness)
%{ 
    This function checks the performance of the objective function around
    the optimal point by plotting slices.

    There is no copyright accompanying this code. Feel free to replicate, 
    post online, or otherwise use as you wish. Please credit the author 
    when you do so. 

    Alexandros Theloudis, LISER & UCL
    Email: a.theloudis@gmail.com

    -----------------------------------------------------------------------
%}

%   Load data:
x0 = phat ;
w0 = what ;

%   Initial statements and switches:
numpoints2plot  = 200 ;             % number of points to plots per parameter
numparams       = length(x0) ;      % number of parameters
mCAvg           = [2.0932;1/2.0932];% pertains to Frisch symmetry        

%   Graphs attributes:
titlex = {'\eta_{h1,w1}','\eta_{h1,w2}','\eta_{h2,w1}','\eta_{h2,w2}'};


%%  1.  INITIALIZE PARAMETER AND FUNCTION VALUE GRID
%   Declare global slice items.
%   -----------------------------------------------------------------------

%   Initial parameter grid:
x0grid = repmat(x0,1,numpoints2plot);

%   Generate matrixes to hold x and f values, and feasibility indicators:
xvalues = zeros(numparams,numpoints2plot);
fvalues = zeros(numparams,numpoints2plot);
feasibl = ones(numparams,numpoints2plot);

%   Parameter bounds:
xb = [  -2.0  2.0 ; ... % eta_h1_w1 
        -2.0  2.0 ; ... % eta_h1_w2 
        -2.0  2.0 ; ... % eta_h2_w1 
        -2.0  2.0]; ... % eta_h2_w2

%   Generate admissible range by parameter:
%   'fineness' is between 0 and 0.99. The highest, the most fine.
if (fineness >= 1) || (fineness < 0)
    error('Input value for finess must be non-negative and strictly less than 1.')
end
xrange = (1-fineness)*0.5*(xb(:,2)-xb(:,1)) ; % fineness is between 0 and 0.99. The highest, the most fine.


%%  2.  LOOP OVER PARAMETERS, CALCULATE FUNCTION, AND PLOT SLICES
%   Per parameter, declare parameter values and loop over. Plot slices.
%   -----------------------------------------------------------------------

for p = 1:1:numparams
    
    %   Generate p-parameter grid keeping everything else constant:
    %clearvars xeval ;
    xeval = x0grid ; 
    xeval(p,:) = xeval(p,:) + [linspace(-xrange(p),0,numpoints2plot/2) linspace(0,xrange(p),numpoints2plot/2)] ;
    
    %   Check grid falls within admissible bounds and if not, update:
    for k = 1:1:numpoints2plot
        if xeval(p,k) < xb(p,1)     % lower bound
            xeval(p,k) = xb(p,1);
        end
        if xeval(p,k) > xb(p,2)     % upper bound
            xeval(p,k) = xb(p,2);
        end
    end
    
    %   Impose linear equality constraints:
    
    %   ! given parameter: eta_h1_w2 (p = 2)
    %   ! eta_h2_w1 = eta_h1_w2 * E[W1*H1/W2*H2]
    if p == 2 
        xeval(p+1,:) = mCAvg(1)*xeval(p,:);
    end
 
    %   ! given parameter: eta_h2_w1 (p = 3)
    %   ! eta_h1_w2 = eta_h2_w1 * E[W2*H2/W1*H1]
    if p == 3
        xeval(p-1,:) = mCAvg(2)*xeval(p,:);
    end

    %   Store x values for the relevant parameter and evaluate function:
    xvalues(p,:) = xeval(p,:); 
    
    %   Evaluate function at each point of the grid:
    for k = 1:1:numpoints2plot
        fvalues(p,k) = laborsupply_structure(xeval(:,k),w0,mD,mInd);
    end

    %   Check that the final point satisfies parameter vector bounds as 
    %   well as non-linear inequality constraints:
    for k = 1:1:numpoints2plot
        
        %   parameter bounds:
        for pp = 1:1:numparams
            if (xeval(pp,k) < xb(pp,1)) || (xeval(pp,k) > xb(pp,2))
                feasibl(p,k) = 0;
            end
        end
    end

    %   Hold on to start plotting:
    hold on
    
    %   Before plotting, check if feasibility is satisfied over domain of 
    %   potential graph. If mean(feasibl) ~= 1, then some x points in the 
    %   parameter graph will not be feasible and must be shaded out:
    
    if mean(feasibl(p,:)) ~= 1
        
        %   Plot starting graph:
        plot(xvalues(p,:),fvalues(p,:),'b');

        %   Obtain existing graph axis limits:
        xLims = get(gca,'XLim');
        yLims = get(gca,'YLim');

        %   Shade entire graph:
        pltpatch = patch([xLims(1) xLims(2) xLims(2) xLims(1)],[yLims(1) yLims(1) (yLims(2)-0.000001) (yLims(2)-0.000001)],'r');
        pltpatch.FaceColor = [0.745098 0.745098 0.745098]; % light grey

        %   Obtain area to un-shade:
        unshaded = feasibl(p,:).*xvalues(p,:); % get non-fesible points in new vector
        unshaded = unshaded(unshaded~=0);      % get rid of zeros (zero = feasible)
        unshaded_xmin = min(unshaded);
        unshaded_xmax = max(unshaded); 

        %   Check if unshaded area is empty and adjust:
        if (isequal(unshaded_xmin,unshaded_xmax) == 1) || (isempty(unshaded_xmin) == 1) || (isempty(unshaded_xmax) == 1)
            unshaded_xmin = phat(p) - (xLims(2)-xLims(1))/200 ;
            unshaded_xmax = phat(p) + (xLims(2)-xLims(1))/200 ;
        end

        %   Unshade:
        patch([unshaded_xmin unshaded_xmax unshaded_xmax unshaded_xmin],[yLims(1) yLims(1) (yLims(2)-0.000001) (yLims(2)-0.000001)],'w');
    end
    
    %   Draw graph (or redraw starting graph as its curve is otherwise 
    %   covered behind shade):
    plt = plot(xvalues(p,:),fvalues(p,:),'b');
    plt.LineWidth = 2;      % curve thickness
    xlabel(titlex(p));      % label x axis
    ylabel('f value');      % label y axis

    %   Plot the estimated parameter value:
    yLims = get(gca,'YLim');
    pmean = plot([phat(p),phat(p)],[yLims(1),(yLims(2)-0.000001)],'r--');
    pmean.LineWidth = 2;  % vertical line thickness

    hold off
    pause
    close all
    
    %   Report progress:
    fprintf('Finished slice of parameter : %u\n',p)
    
end %parameterloop

end