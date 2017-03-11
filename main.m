%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thresholding with cosine edge detection and Hough transform            %
% Using Bilateral filter and rotation                                    %
%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN M-FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function main(filename)
    %     clear all; close all;
    %     filename = 'Bunsen flames new\methane\HIGH RE\High Re phi=1.2\DSC_0346.JPG';
    F = imread(filename); % Read an image
    figure; subplot(1, 3, 1); imshow(F); title('Initial RGB image');

    % Convert rgb image to gray
    B1 = rgb2gray(F);
    subplot(1, 3, 2); imshow(B1); title('Gray image');

    % Crop gray image. Take only flame's ROI
    B = im_cropping(B1);
    subplot(1, 3, 3); imshow(B); title('Cropped image');

    % Calculate max value of cropped image
    Bmax = max(max(double(B)));


    %%%%%%%%%%%%%%%%%%% BILATERAL Filtering %%%%%%%%%%%%%%%%%
    w = 3; % Bilateral filter half-width
    sigma = [3 10]; % Bilateral filter standard deviations
    % Apply bilateral filter to each image.
    B = bfilter2(double(B)/Bmax, w, sigma);
    

    %%%%%%%%%%%%%%%%%%% Find Edges %%%%%%%%%%%%%%%%%
    M = max(B,[],2);

    BB = double(B)./(double(M)*ones(1,size(B, 2)));

    % Find an appropriate threshold value
    thres = thresh_finder(BB);

    CC = zeros(size(BB, 1), size(BB, 2));
    CC(BB >= thres) = 1; % Thresholding
    CC(BB < thres) = 0;

    %%%%%%%%%%%%%%%%%% Extra processing %%%%%%%%%%%%
    % Fill image holes between edges
    CC = imfill(CC, 'holes');
    % Find connected components in binary image
    connected = bwconncomp(CC, 4);
    % Measure area of image regions
    marea = regionprops(connected, 'Area');
    % Remove small objects from binary image
    CC = bwareaopen(CC, max([marea(:).Area]));

    for it = 1:1:11
       CC = medfilt2(CC, [it, it], 'zeros'); 
    end
    
    % Check for close-wise or counterclock-wise rotation
    [l_slope1, r_slope1, slope_sum1] = func_sum_hough(CC, 0);
    fprintf('\n.........................................................................\n');
    fprintf('%s\n', filename);
    fprintf('\nInitially\nleft: %3.2f right: %3.2f sum: %3.2f\n', l_slope1, r_slope1, slope_sum1);
    if slope_sum1 > 180
        l_or_r = 'counter_clock';
    elseif slope_sum1 < 180
        l_or_r = 'clock';
    else
        l_or_r = ' ';
    end


    %%%%%%%%%%%%%%%%%%% ROTATION %%%%%%%%%%%%%%%%%%%
    [m, n] = size(CC);
    % Take a middle region of the image
    rot_CC = CC(m/2-size(B,1)/4:m/2+size(B,1)/4, :);
    stop = 0;
    theta = 0.2; % Initial rotation angle

    % Rotate until sum ~180
    while stop <= 179.5 || stop >= 180.5

        if strcmp(l_or_r, 'counter_clock') == 1
            rot_CC2 = imrotate(rot_CC, theta);
        elseif strcmp(l_or_r, 'clock') == 1
            rot_CC2 = imrotate(rot_CC, -theta);
        else
            fprintf('\nNo need for rotation\n');
            break
        end
        % 
        [l_slope2, r_slope2, slope_sum2] = func_sum_hough(rot_CC2, 0);

        theta = theta + 0.1;
        stop = slope_sum2;

    end

    if strcmp(l_or_r, 'counter_clock') == 1 || strcmp(l_or_r, 'clock') == 1
        fprintf('\n(Rotation) Convergence to 180\nleft: %3.2f right: %3.2f sum: %3.2f\n', l_slope2, r_slope2, slope_sum2);
    end


    PP = CC.*double(B);
    figure; subplot(1, 2, 1); imagesc(PP); title('Detected edges, binary image');


    % Rotate clock or counter-clock
    if strcmp(l_or_r, 'counter_clock') == 1
            CC = imrotate(CC, theta);
    elseif strcmp(l_or_r, 'clock') == 1
            CC = imrotate(CC, -theta);
    end

    subplot(1, 2, 2); imagesc(CC); title('Rotated edges of binary Image');

    % Hough Transform
    [l_slope, r_slope, slope_sum] = func_sum_hough(CC, 1);
    fprintf('\nFinal result (rotated or not)\nleft: %3.2f right: %3.2f sum: %3.2f\n\n', l_slope, r_slope, slope_sum);
    %%%%%%%%%%%%%%%%%%

    % Call curve estimation function
    [step, f_r_top, f_r_mid, f_r_bot, f_l_top, f_l_mid, f_l_bot,...
     df_r_top, df_r_mid, df_r_bot, df_l_top, df_l_mid, df_l_bot, x_r, y_r, x_l, y_l] = fit_intern_edge(CC, l_slope, r_slope);

    % Call function that calculates the area of surface of revolution
    [area_l, area_r] = area_surface_revol(step, f_r_top, f_r_mid, f_r_bot, f_l_top, f_l_mid, f_l_bot,...
                      df_r_top, df_r_mid, df_r_bot, df_l_top, df_l_mid, df_l_bot);

    if(area_l <= area_r)
        fprintf('\nRight and left curve area ratio (left/right): %1.4f\n', area_l/area_r);
    else
        fprintf('\nRight and left curve area ratio (right/left): %1.4f\n', area_r/area_l);
    end

    % Delete the comment so as to let visedge3 function to be executed
%   visedge3(x_r, y_r, x_l, y_l);

    fprintf('*************************************************************************\n');
    fprintf('*************************************************************************\n\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Write results to txt file
% fileID = fopen('text_results.txt','at+');
% filename1 = filename(end-11:end);
% fprintf(fileID, filename1);
% fprintf(fileID, '\n*************************************************************************\n');
% fprintf(fileID, '\nFinal result (rotated or not)\nleft: %3.2f right: %3.2f sum: %3.2f\n\n', l_slope, r_slope, slope_sum);
% fprintf(fileID, 'Left curve area: %d\n', area_l);
% fprintf(fileID, 'Right curve area: %d\n', area_r);
% if(area_l <= area_r)
%     fprintf(fileID, '\nRight and left curve area ratio (left/right): %1.4f\n', area_l/area_r);
% else
%     fprintf(fileID, '\nRight and left curve area ratio (right/left): %1.4f\n', area_r/area_l);
% end
% fprintf(fileID, '*************************************************************************\n');
% fprintf(fileID, '*************************************************************************\n\n');
% fopen('text_results.txt');