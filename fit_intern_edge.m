%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Curve estimation of internal edges and plotting %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input:  im: the image after rotation
%         left_slp: the right slope after rotation (Hough result)
%         right_slp: the left slope after rotation (Hough result)
% Output: step: Distance between points in top and bottom regions
%         f_r_top: f for right top region
%         f_r_mid: f for right middle region
%         f_r_bot: f for right bottom region
%         f_l_top: f for left top region 
%         f_l_mid: f for left middle region
%         f_l_bot: f for left bottom region
%         df_r_top: Derivative of f for right top region
%         df_r_mid: Derivative of f for right middle region
%         df_r_bot: Derivative of f for right bottom region
%         df_l_top: Derivative of f for left top region
%         df_l_mid: Derivative of f for left middle region
%         df_l_bot: Derivative of f for left middle region
%         x_r: Right x points
%         y_r: Right y points 
%         x_l: Left x points
%         y_l: Left y points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [step, f_r_top, f_r_mid, f_r_bot, f_l_top, f_l_mid, f_l_bot,...
         df_r_top, df_r_mid, df_r_bot, df_l_top, df_l_mid, df_l_bot, x_r, y_r, x_l, y_l] = fit_intern_edge(im, left_slp, right_slp)

    warning('off','all')
    warning

%     im = CC; % input  , after rotation image 
%     left_slp = l_slope; % input,                
%     right_slp = r_slope; % input,               
    
    step = 3; % Point selection step
    
    % Calculate the limits of x between which there are 4 transitions 1<->0
    dfs1 = zeros(size(im, 1)-1, size(im, 2)-1);
    for i = 1:size(im, 1)
        dfs1(i, :) = diff(im(i, :));
    end
    
    num_nnz = zeros(size(dfs1, 1), 1);
    for i = 1:size(num_nnz, 1)
        num_nnz(i) = nnz(dfs1(i, :));
    end
    
    % Limits of x between which there are 4 transitions 1<->0 
    top_lim = max(find(num_nnz == 4, 1, 'first'));
    bottom_lim =  min(find(num_nnz == 4, 1, 'last'));
    

    % Initialize coordinate accumulator vectors for the right line
    x_r = zeros(1, ceil((bottom_lim - top_lim)/step)); 
    y_r = zeros(1, ceil((bottom_lim - top_lim)/step));
    % Initialize coordinate accumulator vectors for the left line
    x_l = zeros(1, ceil((bottom_lim - top_lim)/step));
    y_l = zeros(1, ceil((bottom_lim - top_lim)/step));

    pointer = 1;
    
   
   % Collect the points and store their coordinate values 
   % Search between top_lim and bottom_lim
   for iter = top_lim:step:bottom_lim % 
        
        dfs = diff(im(iter, :)); % Find differences
        % figure; stem(dfs);
    
        % Left line point detection
        tmp1 = find(dfs == -1, 1, 'first'); % The first -1 of dfs
        
        % If a point is found store values of its coordinates
        if (isempty(tmp1) == 0)
            x_l(1, pointer) = iter;
            y_l(1, pointer) = tmp1;
        end

        % Right line point detection
        tmp2 = find(dfs == 1, 1, 'last');  % The last 1 of dfs

        % If a point is found store values of its coordinates
         if (isempty(tmp2) == 0)
            x_r(1, pointer) = iter;
            y_r(1, pointer) = tmp2;
         end

        pointer = pointer + 1;
        
   end
   
   % Find those x with angle close to the angles (left or right) returned by Hough 
   Z_l = atan(diff(x_l)./diff(y_l))*180/pi;
   Z_r = atan(diff(x_r)./diff(y_r))*180/pi;
   
   Z_l = abs(Z_l);
   Z_r = 180 - Z_r; 
   
   
   [counts1, centers1] = hist(Z_l, 0 : 0.1 : 179-0.1);
   % figure; bar(centers1,counts1);
   [counts2, centers2] = hist(Z_r, 0 : 0.1 : 179-0.1);
   % figure; bar(centers2,counts2);
   
   clear ind_l_mid; clear ind_r_mid;    % 
   
   % [~, temp1] = max(counts1(counts1 ~= max(counts1)));
   [~, temp1] = max(counts1);
   ind_l_mid = find(Z_l < temp1*0.1 + 1 & Z_l > temp1*0.1 - 1);
   % [~, temp2] = max(counts2(counts2 ~= max(counts2)));
   [~, temp2] = max(counts2);
   ind_r_mid = find(Z_r < temp2*0.1 + 1 & Z_r > temp2*0.1 - 1);
   
   
   % Sample middle points to avoid very small top and bottom regions
   ind_l_mid = ind_l_mid(2:1:end-1);
   ind_r_mid = ind_r_mid(2:1:end-1);

   
   % Indexes of bottom zone points
   ind_l_bot = find(x_l>=x_l(ind_l_mid(end)));
   ind_r_bot = find(x_r>=x_r(ind_r_mid(end)));
   
   % Indexes of top zone points
   ind_l_top = find(x_l<=x_l(ind_l_mid(1)));
   ind_r_top = find(x_r<=x_r(ind_r_mid(1)));
   
   figure; imshow(im); hold on;
   plot(y_l(ind_l_top), x_l(ind_l_top),'*r','MarkerSize',3); hold on;
   plot(y_r(ind_r_top), x_r(ind_r_top),'*r','MarkerSize',3); hold on;
   plot(y_l(ind_l_bot), x_l(ind_l_bot),'*g','MarkerSize',3); hold on;
   plot(y_r(ind_r_bot), x_r(ind_r_bot),'*g','MarkerSize',3); hold on;
   plot(y_l(ind_l_mid), x_l(ind_l_mid),'*b','MarkerSize',3); hold on;
   plot(y_r(ind_r_mid), x_r(ind_r_mid),'*b','MarkerSize',3);
   axis ij; % Change origin
   
   % Shift y values in order to be located around 0
   tmp_peak = (y_r(1) + y_l(1))/2;
   y_r = y_r - tmp_peak;
   y_l = y_l - tmp_peak;
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%% PIECEWISE CURVE FITTING & PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   figure;
   plot(y_l, x_l, '+r', y_r, x_r, '+r',... 
        y_l(ind_l_mid), x_l(ind_l_mid), 'og', y_r(ind_r_mid), x_r(ind_r_mid), 'og'); hold on;
   axis ij; % Change origin 
   fprintf('\nPiecewise curve fitting');
   
    % First call curve_constr_fit() function for middle regions
    % calculating constr. points in derivative. These points are being
    % propagated to the calls of the function for bot and top regions as
    % fixed constr. points
   
    % Call curve contraint fit for left middle part
    % 3rd input is set to 0. No contr. point to be propagated
    [f_l_mid, df_l_mid, deg, new_ind_l_mid] = curve_constr_fit(x_l, y_l, ind_l_mid, 0, 'mid');
    fprintf('\n\nLeft middle polynomial degree: %d\n', deg);
    
    plot(f_l_mid, x_l(new_ind_l_mid), '-b'); hold on;
    axis ij
    
    % Call curve contraint fit for right middle part
    % 3rd input is set to 0. No contr. point to be propagated
    [f_r_mid, df_r_mid, deg, new_ind_r_mid] = curve_constr_fit(x_r, y_r, ind_r_mid, 0, 'mid');
    fprintf('Right middle polynomial degree: %d\n', deg);
    
    plot(f_r_mid, x_r(new_ind_r_mid), '-b'); hold on;
    axis ij
    
    
    % Call curve contraint fit for left bottom part
    [f_l_bot, df_l_bot, deg] = curve_constr_fit(x_l, y_l, ind_l_bot, f_l_mid(end),'bot');
    fprintf('Left bottom polynomial degree: %d\n', deg);
    
    plot(f_l_bot, x_l(ind_l_bot), '-b'); hold on;
    axis ij
    
    % Call curve contraint fit for right bottom part
    [f_r_bot, df_r_bot, deg] = curve_constr_fit(x_r, y_r, ind_r_bot, f_r_mid(end), 'bot');
    fprintf('Right bottom polynomial degree: %d\n', deg);
    
    plot(f_r_bot, x_r(ind_r_bot), '-b'); hold on;
    axis ij
        
    % Call curve contraint fit for left top part
    [f_l_top, df_l_top, deg] = curve_constr_fit(x_l, y_l, ind_l_top, f_l_mid(1), 'top');
    fprintf('Left top polynomial degree: %d\n', deg);
      
    % Call curve contraint fit for right top part
    [f_r_top, df_r_top, deg] = curve_constr_fit(x_r, y_r, ind_r_top, f_r_mid(1), 'top');
    fprintf('Right top polynomial degree: %d\n', deg);
    
    % Plot left top part and right top part simultaneously
    plot([f_l_top(end:-1:1)' f_r_top'], [x_l(ind_l_top(end:-1:1)) x_r(ind_r_top)], '-b'); hold on;
    
    peak_y = (f_l_top(1) + f_r_top(1))/2;
    plot(peak_y, x_l(1), '*y');
    axis ij   
    
    title('Points of internal edges'); legend('All points', '', 'Hough slope points', '', 'Constraint points');   
end    