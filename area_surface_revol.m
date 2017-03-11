%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the area of surface of revolution using left and right  %
% fitted edge curves                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input: step: Distance between points in top and bottom regions
%        f_r_top: f for right top region
%        f_r_mid: f for right middle region
%        f_r_bot: f for right bottom region
%        f_l_top: f for left top region 
%        f_l_mid: f for left middle region
%        f_l_bot: f for left bottom region
%        df_r_top: Derivative of f for right top region
%        df_r_mid: Derivative of f for right middle region
%        df_r_bot: Derivative of f for right bottom region
%        df_l_top: Derivative of f for left top region
%        df_l_mid: Derivative of f for left middle region
%        df_l_bot: Derivative of f for left middle region
% Output: sum_l: Left curve area
%         sum_r: Right curve area
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sum_l, sum_r] = area_surface_revol(step, f_r_top, f_r_mid, f_r_bot, f_l_top, f_l_mid, f_l_bot,...
                            df_r_top, df_r_mid, df_r_bot, df_l_top, df_l_mid, df_l_bot)
                        
    % Piecewise area calculation for left and right curves
    
    % Init three region counters
    ss1 = 0;  ss2 = 0; ss3 = 0;
            
    % Calcultion of area for right top region
    for i = 2:size(f_r_top, 1)
        % dx is equal to step
        ss1 = ss1 + 2.*pi.*f_r_top(i).*sqrt(1 + df_r_top(i).^2) * step; 
    end
    
    % Calcultion of area for right mid region
    for i = 2:size(f_r_mid, 1)
        % 
        ss2 = ss2 + 2.*pi.*f_r_mid(i).*sqrt(1 + df_r_mid(i).^2) * step;% Dx(end-i+2); 
    end

    % Calcultion of area for right bot region
    for i = 2:size(f_r_bot, 1)
        % dx is equal to step
        ss3 = ss3 + 2.*pi.*f_r_bot(i).*sqrt(1 + df_r_bot(i).^2) * step; 
    end

    % Summary of area of the right curve  
    sum_r = ss1 + ss2 + ss3;
    
    % Clear counters
    ss1 = 0; ss2 = 0; ss3 = 0;
    clear Dx

    % Calcultion of area for left top region
    for i = 2:size(f_l_top, 1)
        % dx is equal to step
        ss1 = ss1 + abs(2.*pi.*f_l_top(i).*sqrt(1 + df_l_top(i).^2)) * step; 
    end

    % Calcultion of area for left middle region
    for i = 2:size(f_l_mid, 1)
        % 
        ss2 = ss2 + abs(2.*pi.*f_l_mid(i).*sqrt(1 + df_l_mid(i).^2)) * step; % Dx(end-i+2); 
    end

    % Calcultion of area for left bot region
    for i = 2:size(f_l_bot, 1)
        % dx is equal to step
        ss3 = ss3 + abs(2.*pi.*f_l_bot(i).*sqrt(1 + df_l_bot(i).^2)) * step; 
    end

    % Summary of area of the left curve
    sum_l = ss1 + ss2 + ss3;
    fprintf('\n\nLeft curve area: %d\n', sum_l);
    fprintf('Right curve area: %d\n', sum_r);
    
end