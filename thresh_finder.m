%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the ideal threshold value %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input: BB: Image for processing
% Output: ret: Calculated threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ret = thresh_finder(BB)

    backup_BB = BB; % Backup of BB

    % Initialize theta step
    th_step = 0.25;
    % Thresholds range
    all_thres = 0.95:-0.01:0.85;
    % Line number and length vectors
    num_lns = zeros(1, size(all_thres, 2));
    lng_lns = zeros(1, size(all_thres, 2));
    i = 0; % Iter
    % Initial threshold
    thres = 0.95;

    % Find edges for thresholds in range [0.85, 0.95]
    % Select this threshold providing edges without discontinueties
    % as thick as posible
    while thres > 0.84

        BB = (backup_BB>thres).*backup_BB; 
        CC = (cos(pi*BB));
        CC = abs((CC<0).*CC);

        % Calculate gradient of image
        [x y] = gradient(double(CC));
        % Apply Hough Transform
        [H2, theta2, rho2] = hough(abs(x), 'Theta', -90:th_step:(90-th_step)); 
        Peaks  = houghpeaks(H2, 6, 'threshold', ceil(0.9*max(H2(:))));
        lines = houghlines(CC, theta2, rho2, Peaks, 'FillGap', 5, 'MinLength', 10);
        % figure, imshow(CC), hold on

        max_len = 0;
        for k = 1:length(lines)
            xy = [lines(k).point1; lines(k).point2];
            
            % Determine the endpoints of the longest line segment
            len = norm(lines(k).point1 - lines(k).point2);
            if(len > max_len)
                max_len = len;
            end
        end

        % highlight the longest line segment
        % plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','blue');


        i = i + 1;
        thres = thres - 0.01; % Decrease threshold
        num_lns(i) = size(lines, 2); % Store # of lines
        lng_lns(i) = max_len; % Store length of biggest line
        
    end
    
    % Find threshold with min formulated lines. From those choose the one
    %  having the longest line
    % num_lns
    % lng_lns
    index_1 = find(num_lns == min(num_lns));
    [i, index_2] = max(lng_lns(index_1));
    
    ret = all_thres(index_1(index_2));
end