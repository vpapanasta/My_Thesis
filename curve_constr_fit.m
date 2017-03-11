%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit curves with constraints on derivative using lsqlin() function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input: x_c: x points of left or right part
%        y_c: y points of left or right part
%        ind: indeces of the specific region
%        y_pr_constr: propagated y constraint point of other region(0-value if no need for propagation)
%        type: the type of the region (mid, bot, top)
% Output: ret1: contains the final fitted points
%         ret2: contains the derivative of the final fitted points
%         ret3: return polynomial degree
%         ret4: new middle region index (for all mid points)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ret1, ret2, ret3, ret4] = curve_constr_fit(x_c, y_c, ind, y_pr_constr, type)
    
    % Use all points between ind(1) and ind(end) 
    x = x_c(find(x_c >= x_c(ind(1)) & x_c <= x_c(ind(end))));
    y = y_c(find(x_c >= x_c(ind(1)) & x_c <= x_c(ind(end))));
    
    % Constraint Points
    x0 = x_c(ind(end));
    x1 =  x_c(ind(1));
    y0 = y_c(find(x_c == x_c(ind(end))));
    y1 = y_c(find(x_c == x_c(ind(1))));

    x = x(:); % Reshape the data into a column vector
    y = y(:);
   
    
    % Estimate degree 'n' of polynomial to fit using Bayesian information criterion (BIC) 
    n = degree_bic(x_c, y_c, ind);
        
    %%%%%%%%%%%%%%%%%%%%%%%% Initial p-values calculation %%%%%%%%%%%%%%%%%%%%%%%%
    
    % Construct Vandermode matrix V for x   
    V(:,n) = ones(length(x),1,class(x));
    for j = n-1:-1:1
         V(:,j) = x.*V(:,j+1);
    end

    % Use linear equality constraints to force the curve to hit the required point. In
    % this case, 'Aeq(1, :)' is the Vandermoonde matrix for 'x0' and
    % 'Aeq(2, :)' for 'x1'. 
    Aeq = zeros(2, n);
    Aeq(1, :) = x0.^(n-1:-1:0);
    Aeq(2, :) = x1.^(n-1:-1:0);

    % and 'beq' is the value the curve should take at that point
    beq = zeros(1, 2);
    beq(1) = y0;
    beq(2) = y1;
    

    % Solve constrained linear least-squares problem calculating parameters p
    options = optimset('algorithm','interior-point'); 
    options.Display = 'off';
    %
    p = lsqlin(V, y, [], [], Aeq, beq, [], [], [], options);

    
    
   %%%%%%%%%%%%%%%%%%%%%%%%%% End of initial p-values calculation%%%%%%%%%%%%%%%%%%%    
    
    % Make Constraints on the derivatives
    % If middle region calculate A1, A2 on derivative. If bottom region not
    % derivative for both. If top region derivative for A1 only.
    if isequal(type, 'mid') % If middle region
        A1 = [(n-2:-1:1).*x0.^(n-3:-1:0) 0 0]; 
        A2 = [(n-2:-1:1).*x1.^(n-3:-1:0) 0 0];
        Aeq = [A1; A2];
    elseif isequal(type, 'bot') % If bottom region
        % I need constraint only for the connection with mid region. End
        % point is not constrained
        %A1= x0.^(n-1:-1:0);
        A2= x1.^(n-1:-1:0);
        Aeq = A2; 
    else % If top region
        A1= x0.^(n-1:-1:0);
        A2 = [(n-2:-1:1).*x1.^(n-3:-1:0) 0 0]; 
        Aeq = [A1; A2];
    end
    

    % Calculate derivatives der_y0, der_y1 of x0, x1
    a = p;
    if isequal(type, 'mid') % If middle region
        V_x0 = A1;
        V_x1 = A2;
    elseif isequal(type, 'top') % If top region
        V_x1 = A2;
    end
    
    % If the region is a middle one calculate der_y0 and der_y1. Else if bottom region, der_y0 equals
    % to middle region's der_y1 (propagate) and for the down point we use
    % y1 value. Else if top region, der_y1 equals to middle region's der_y0 and for the up point der_y0 is calculated 
    if isequal(type, 'mid') % If middle region
        der_y0 = V_x0 * a;
        der_y1 = V_x1 * a;
        beq = [der_y0, der_y1]; %
    elseif isequal(type, 'bot') % If bottom region
        % I need constraint only for the connection with mid region.
        % beq = [y0, y_pr_constr]; %
        beq = y_pr_constr;
    else % If top region
        der_y1 = V_x1 * a;
        beq = [y_pr_constr, der_y1]; %
    end
    

    % Call LSQLIN with options to prevent warnings
    opts = optimset('algorithm', 'interior-point'); 
    opts.LargeScale = 'off'; 
    opts.Display = 'off';
    p = lsqlin(V, y ,[], [], Aeq, beq, [], [], [], opts);

    % Use POLYVAL to evaluate the fitted curve (f and derivative of f)
    f = polyval(p, x);
    der_f = polyval(polyder(p), x);
    
    
    % Plot point to go through
    plot(y0,x0,'k*','linewidth', 1);
    plot(y1,x1,'k*','linewidth', 1); 
    
    % If the region is a mid or top one plot its derivative contr point with cyan color
    if isequal(type, 'mid')
        plot(f(end),x0,'c*','linewidth', 1);
        plot(f(1), x1,'c*','linewidth', 1);
    elseif isequal(type, 'top')
        plot(f(1), x1,'c*','linewidth', 1);
    end
    
    ret1 = f; % Return the fitted points
    ret2 = der_f;
    ret3 = n; % Return polynomial degree
    ret4 = find(x_c >= x_c(ind(1)) & x_c <= x_c(ind(end)));
    
end