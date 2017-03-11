%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Flame's slopes calculation using Hough Transform %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input: image: Image for processing
%        show: Show flag
% Output: ret1: Left slope
%         ret2: Right slope
%         ret3: Summary of the two slopes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ret1, ret2, ret3] = func_sum_hough(image, show)
    % Initialize theta step
    th_step = 0.1;
    
    % Calculate gradient of image
    [x y] = gradient(double(image));
    
    % Apply Hough Transform
    [H2, theta2, rho2] = hough(abs(x), 'Theta', -90:th_step:(90-th_step)); 
    
    if show == 1
       figure; imagesc(H2); title('x gradient -> Hough');
       figure; mesh(H2); title('x gradient -> Hough');
       figure; imshow(imadjust(mat2gray(H2))); title('Hough Transform of image');
       colormap(hot);
    end

    % Find max angles with max accumulated values (left & right)    
    [C, l_slope] = max(max(H2(:,1:90/th_step))); 
    [C, r_slope] = max(max(H2(:,(90/th_step+1):end)));
    
    % Slopes' correction
    l_slope = l_slope*th_step;
    r_slope = r_slope*th_step;
    
    r_slope = r_slope + 90;
    % Slopes' summary
    slope_sum = r_slope + l_slope;
    
    % Return values
    ret1 = l_slope;
    ret2 = r_slope;
    ret3 = slope_sum; 

end