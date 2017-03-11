%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation of optimal polynomial degree using Bayesian Information %
% Criterion                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input: x_c: x points of left or right part
%        y_c: y points of left or right part
%        ind: indeces of the specific region
% Output: ret: the degree of polynomial for which BIC is minimized 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ret = degree_bic(x_c, y_c, ind)

    if (ind ~= 999) % In case of piecewise fitting
        % The sample points
        x = x_c(ind);
        y = y_c(ind);

        % Constraint Points
        x0 = x_c(ind(end));
        x1 = x_c(ind(1));
        y0 = y_c(find(x_c == x_c(ind(end))));
        y1 = y_c(find(x_c == x_c(ind(1))));

        x = x(:); % Reshape the data into a column vector
        y = y(:);
    else % Non piecewise fitting
        % Constraint Points
        x0 = x_c(end);
        x1 = x_c(1);
        y0 = y_c(end);
        y1 = y_c(1);

        x = x_c(:);
        y = y_c(:);
    end
    
    test_num = 50; % # of degree tests
    n = size(x, 1); % # of observations
    bic = zeros(test_num, 1); % Init bic vector

    % For each degree
    for deg = 1:test_num
        
        % Construct Vandermode matrix V for x   
        V(:,deg+1) = ones(length(x),1,class(x));
        for j = deg:-1:1
             V(:,j) = x.*V(:,j+1);
        end

        % Use linear equality constraints to force the curve to hit the required point. In
        % this case, 'Aeq(1, :)' is the Vandermoonde matrix for 'x0' and
        % 'Aeq(2, :)' for 'x1'.
        Aeq = zeros(2, deg+1);
        Aeq(1, :) = x0.^(deg:-1:0);
        Aeq(2, :) = x1.^(deg:-1:0);

        % and 'beq' is the value the curve should take at that point
        beq = zeros(1, 2);
        beq(1) = y0;
        beq(2) = y1;


        % Solve constrained linear least-squares problem calculating parameters p
        options = optimset('algorithm','interior-point'); % (rec)
        % options = optimset('algorithm','trust-region-reflective'); options.LargeScale = 'on';
        % options = optimoptions(@lsqlin,'Algorithm','active-set');  options.LargeScale = 'off';
        options.Display = 'off';
        [p, ~, Residual, ~, ~, ~] = lsqlin(V, y, [], [], Aeq, beq, [], [], [], options);

        % Error equal to residual from lsqlin() function 
        Error = Residual;
        
        % Sum Squared Error calculation
        sse = sum(Error.^2);
        
        % Bayesian information criterion calculation        
        bic(deg) = n*log(sse/n) + deg*log(n);

    end
    
    %     figure; plot(1:test_num, bic, '-xb');
    %     xlabel('Degree');
    %     ylabel('BIC');
    
    % Return this degree 'deg' for which BIC is minimized
    [~,ret] = min(bic(find(bic~=-inf))); % Don't take into consideration -inf values
    if ret == 1
        ret = 2;
    end
    
end