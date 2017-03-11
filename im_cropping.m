%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Flame Image Cropping %
%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input: grayImage: Image for processing
% Output: ret: Cropped image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ret = im_cropping(grayImage)
    % im must be a gray image

    % Calculate Global image threshold using Otsu's method
    level = graythresh(grayImage);
    % Converts gray image to binary
    binaryImage = im2bw(grayImage, level);
    % Filter binary image in order to remove salt&pepper noise. Use of 7x7
    % window for better results
    binaryImage = medfilt2(binaryImage, [7 7]);

    % Find horizontal and vertical point to specify flame's ROI
    horizontal = any(binaryImage, 1);
    vertical = any(binaryImage, 2);

    column1 = find(horizontal, 1, 'first');
    column2 = find(horizontal, 1, 'last');
    row1 = find(vertical, 1, 'first');
    row2 = find(vertical, 1, 'last');

    % Crop image 
    croppedImage = grayImage(row1:row2, column1:column2);
    % Return cropped image
    ret = croppedImage;

end