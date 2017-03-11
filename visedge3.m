%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cone of flame visualization rotating the the right and left fitted %
% edges                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input: x_r: x points of right part
%        y_r: y points of right part
%        x_l: x points of left part
%        y_l: y points of left part
% Output: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function visedge3(x_r, y_r, x_l, y_l)

    gds = 1;
    
    x = abs(x_r);
    y = abs(y_r);
    
    data_r = zeros(size(x(:), 1), 3, 361);
    data_l = zeros(size(x(:), 1), 3, 361);
    figure;
    
    while gds < 3

        %Vertices matrix
        V = [x(:) y(:) zeros(size(y(:)))];
        V_centre = mean(V,1); % Centre, of line
        Vc = V-ones(size(V,1),1)*2*V_centre; % Centering coordinates

        step = 1;

        a = 1; % Angle in degrees

        for i = 0:360/step
            
            i = i + 1;
            a_rad = ((a*pi)./180); % Angle in radians
            E = [0  0 a_rad]; % Euler angles for X,Y,Z-axis rotations

            R = [cos(E(3))  0        sin(E(3));...
                   0        1        0;...
                 -sin(E(3)) 0        cos(E(3))]; %Y-axis rotation
       
            Vrc = [R*Vc']'; % Rotating centred coordinates

            Vr = Vrc+ones(size(V,1),1)*2*V_centre; % Centering coordinates [option 1]
            
            if gds == 1
                data_r(:, 1, i) = Vr(:,1);
                data_r(:, 2, i) = Vr(:,2);
                data_r(:, 3, i) = Vr(:,3);
                plot_r = plot3(data_r(:, 1, i), data_r(:, 2, i), data_r(:, 3, i),'r-', 'MarkerSize', 4); % Rotated around centre of line (right)
            elseif gds == 2
                data_l(:, 1, i) = Vr(:,1);
                data_l(:, 2, i) = Vr(:,2);
                data_l(:, 3, i) = Vr(:,3);
                plot_l = plot3(data_l(:, 1, i), data_l(:, 2, i), data_l(:, 3, i),'b-', 'MarkerSize', 4); % Rotated around centre of line (left)
            end

            xlabel('x'); ylabel('y'); zlabel('z');
            hold all 

            a = a + step;
        %     
        end
        
         gds = gds + 1;
         x = abs(x_l);
         y = abs(y_l);
    end
    legend([plot_r, plot_l], 'Rotated Right Curve', 'Rotated Left Curve');
    
    newd1_r = reshape(data_r(:, 1, :), [size(data_r, 1), size(data_r, 3)]);
    newd2_r = reshape(data_r(:, 2, :), [size(data_r, 1), size(data_r, 3)]);
    newd3_r = reshape(data_r(:, 3, :), [size(data_r, 1), size(data_r, 3)]);
    
    newd1_l = reshape(data_l(:, 1, :), [size(data_l, 1), size(data_l, 3)]);
    newd2_l = reshape(data_l(:, 2, :), [size(data_l, 1), size(data_l, 3)]);
    newd3_l = reshape(data_l(:, 3, :), [size(data_l, 1), size(data_l, 3)]);
    
    figure; contour3(newd1_r, newd3_r, newd2_r, 360); title('Contour3 Right');
    figure; contour3(newd1_l, newd3_l, newd2_l, 360); title('Contour3 Left');
    
    figure; surf(newd1_r, newd3_r, newd2_r); shading interp; title('Surf Right');
    figure; surf(newd1_l, newd3_l, newd2_l); shading flat; title('Surf Left');
    
end