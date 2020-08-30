%% Canny Edge Detection
% Bilal Ahmad

addpath /Users/jeborah/Desktop/

%% Produce document figures
subplot (1,3,1)
lanFlowers = imread("Flowers.jpg");
imshow(lanFlowers);
title('Original Image');
subplot(1,3,2)
CannyEdgeDet("Flowers.jpg", 1.5, .075, .175);
title('\sigma = 1.5, \tau_l = .075, \tau_h = .175');
subplot(1,3,3)
CannyEdgeDet("Flowers.jpg", 3, .075, .175);
title('\sigma = 3.0, \tau_l = .075, \tau_h = .175');
saveas(gcf, 'Flowers', 'fig')

subplot (1,3,1)
lanFlowers = imread("Syracuse_01.jpg");
imshow(lanFlowers);
title('Original Image');
subplot(1,3,2)
CannyEdgeDet("Syracuse_01.jpg", 1.5, .075, .175);
title('\sigma = 1.5, \tau_l = .075, \tau_h = .175');
subplot(1,3,3)
CannyEdgeDet("Syracuse_01.jpg", 1.5, .075, .375);
title('\sigma = 1.5, \tau_l = .075, \tau_h = .375');
saveas(gcf, 'Syr1', 'fig')

subplot (1,3,1)
lanFlowers = imread("Syracuse_02.jpg");
imshow(lanFlowers);
title('Original Image');
subplot(1,3,2)
CannyEdgeDet("Syracuse_02.jpg", 1.5, .125, .175);
title('\sigma = 1.5, \tau_l = .125, \tau_h = .175');
subplot(1,3,3)
CannyEdgeDet("Syracuse_02.jpg", 1.5, .005, .175);
title('\sigma = 1.5, \tau_l = .005, \tau_h = .175');
saveas(gcf, 'Syr2', 'fig')

subplot (1,2,1)
lanFlowers = imread("cctv.jpg");
imshow(lanFlowers);
title('Original Image');
subplot(1,2,2)
CannyEdgeDet("cctv.jpg", .5, .075, .175);
title('\sigma = .5, \tau_l = .075, \tau_h = .175');
saveas(gcf, 'cctv', 'fig')
%% Function Reference - CannyEdgeDet
% This function applies the Canny Edge Detection algorithm to an image
%
% Inputs:
%
% * asImage - name of the image file to perform algorithm on
% * arSigma - sigma value used for initial Gaussian filtering
% * arLowThresh - value of the lower threshold limit
% * arHighThresh - value of upper threshold limit
function CannyEdgeDet(asImage, arSigma, arLowThresh, arHighThresh)
% Import Image
addpath /Users/jeborah/Desktop/cis/691/hw2/code
lanImage = imread(asImage);

% convert image to greyscale if necessary
if (size(lanImage, 3) ~= 1)
    lanImage = rgb2gray(lanImage);
end
%% Calculate Strength and Orientation Images
lrSigma = arSigma;

% imgaussfilt performs a convolution with a 5x5 Gaussian kernel
larSmoothedImage = imgaussfilt(lanImage, lrSigma);

% Compute the horizontal and vertical gradients
lanGradX = [-1, 0, 1; -2, 0, 2; -1, 0, 1];
lanGradY = [1, 2, 1; 0, 0, 0; -1, -2, -1];
larGradImgX = conv2(larSmoothedImage, lanGradX, 'same');
larGradImgY = conv2(larSmoothedImage, lanGradY, 'same');

% Get the orientation image, convert radians to degrees
larOrImage = atan2 (larGradImgY, larGradImgX);
larOrImage = larOrImage*180/pi;

% Get the strength image
larStrImage = (larGradImgX.^2) + (larGradImgY.^2);
larStrImage = sqrt(larStrImage);

%% Apply Non-max Suppression
lnNumRows = size(larSmoothedImage,1);
lnNumCols = size(larSmoothedImage,2);

% make sure all the values in the orientation image are positive before
% beginning non-max suppression
for i=1:lnNumRows
    for j=1:lnNumCols
        if (larOrImage(i,j) < 0)
            larOrImage(i,j) = larOrImage(i,j) + 360;
        end
    end
end

% Update the values in the orientation image to be 0, 45, 90, or 135
for i = 1  : lnNumRows
    for j = 1 : lnNumCols
        if ((larOrImage(i, j) >= 0 ) && ...
                (larOrImage(i, j) < 22.5) || ...
             (larOrImage(i, j) >= 157.5) && ...
                (larOrImage(i, j) < 202.5) || ...
             (larOrImage(i, j) >= 337.5) && ...
                (larOrImage(i, j) <= 360))
            larOrImage(i, j) = 0;
        elseif ((larOrImage(i, j) >= 22.5) && ...
                    (larOrImage(i, j) < 67.5) || ...
                (larOrImage(i, j) >= 202.5) && ...
                    (larOrImage(i, j) < 247.5))
            larOrImage(i, j) = 45;
        elseif ((larOrImage(i, j) >= 67.5 && ...
                    larOrImage(i, j) < 112.5) || ...
                (larOrImage(i, j) >= 247.5 && ...
                    larOrImage(i, j) < 292.5))
            larOrImage(i, j) = 90;
        elseif ((larOrImage(i, j) >= 112.5 && ...
                    larOrImage(i, j) < 157.5) || ...
                (larOrImage(i, j) >= 292.5 && ...
                    larOrImage(i, j) < 337.5))
            larOrImage(i, j) = 135;
        end
    end
end

% Perform the non-maximum supression
larNonMaxSup = zeros (lnNumRows, lnNumCols);

for i=2:lnNumRows-1
    for j=2:lnNumCols-1
        if (larOrImage(i,j)==0)
            larNonMaxSup(i,j) = ...
                (larStrImage(i,j) == max([larStrImage(i,j), ...
                larStrImage(i,j+1), ...
                larStrImage(i,j-1)]));
        elseif (larOrImage(i,j)==45)
            larNonMaxSup(i,j) = ...
                (larStrImage(i,j) == max([larStrImage(i,j), ...
                larStrImage(i+1,j-1), ...
                larStrImage(i-1,j+1)]));
        elseif (larOrImage(i,j)==90)
            larNonMaxSup(i,j) = ...
                (larStrImage(i,j) == max([larStrImage(i,j), ...
                larStrImage(i+1,j), ...
                larStrImage(i-1,j)]));
        elseif (larOrImage(i,j)==135)
            larNonMaxSup(i,j) = ...
                (larStrImage(i,j) == max([larStrImage(i,j), ...
                larStrImage(i+1,j+1), ...
                larStrImage(i-1,j-1)]));
        end
    end
end

larNonMaxSup = larNonMaxSup.*larStrImage;

%% Perform Hysteresis Thresholding

% upper and lower threshold values
lrThreshLow = arLowThresh;
lrThreshHigh = arHighThresh;
lrThreshLow = lrThreshLow * max(max(larNonMaxSup));
lrThreshHigh = lrThreshHigh * max(max(larNonMaxSup));

larHysThresh = zeros(lnNumRows, lnNumCols);

for i = 1  : lnNumRows
    for j = 1 : lnNumCols
        if (larNonMaxSup(i, j) < lrThreshLow)
            larHysThresh(i, j) = 0;
        elseif (larNonMaxSup(i, j) > lrThreshHigh)
            larHysThresh(i, j) = 1;
        elseif (larNonMaxSup(i+1,j)>lrThreshHigh || ...
                larNonMaxSup(i-1,j)>lrThreshHigh || ...
                larNonMaxSup(i,j+1)>lrThreshHigh || ...
                larNonMaxSup(i,j-1)>lrThreshHigh || ...
                larNonMaxSup(i-1, j-1)>lrThreshHigh || ...
                larNonMaxSup(i-1, j+1)>lrThreshHigh || ...
                larNonMaxSup(i+1, j+1)>lrThreshHigh || ...
                larNonMaxSup(i+1, j-1)>lrThreshHigh)
            larHysThresh(i,j) = 1;
        end
    end
end

% cast the values into integers
larResult = 255 - uint8(larHysThresh.*255);
imshow(larResult, []);
end

