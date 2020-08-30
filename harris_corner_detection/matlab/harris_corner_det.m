%% Harris Corner Detection
% Bilal Ahmad 2/21/20

addpath /Users/jeborah/Desktop/
%% Produce Document Figures
subplot (1,3,1)
lanCheckerboard = imread("CheckerBoard.jpg");
imshow(lanCheckerboard);
title('Original Image');
subplot(1,3,2)
CornerDetection("CheckerBoard.jpg", 1.5, 1, 7500);
title('\sigma = 1.5, Q = 3x3, \tau = 7500');
subplot(1,3,3)
CornerDetection("CheckerBoard.jpg", 3, 1, 7500);
title('\sigma = 3, Q = 3x3, \tau = 7500');
saveas(gcf, 'CheckSig', 'fig')

subplot (1,3,1)
lanCheckerboard = imread("CheckerBoard.jpg");
imshow(lanCheckerboard);
title('Original Image');
subplot(1,3,2)
CornerDetection("CheckerBoard.jpg", 1.5, 1, 7500);
title('\sigma = 1.5, Q = 3x3, \tau = 7500');
subplot(1,3,3)
CornerDetection("CheckerBoard.jpg", 1.5, 3, 5000);
title('\sigma = 1.5, Q = 7x7, \tau = 7500');
saveas(gcf, 'CheckNeighborhood', 'fig')

subplot (1,3,1)
lanCheckerboard = imread("Building1.jpg");
imshow(lanCheckerboard);
title('Original Image');
subplot(1,3,2)
CornerDetection("Building1.jpg", 1.5, 1, 7500);
title('\sigma = 1.5, Q = 3x3, \tau = 7500');
subplot(1,3,3)
CornerDetection("Building1.jpg", 1.5, 1, 10000);
title('\sigma = 1.5, Q = 3x3, \tau = 10000');
saveas(gcf, 'BuildingThresh', 'fig')

%% Function Reference - CornerDetection
% This function applies the Canny Edge Detection algorithm to an image
%
% Inputs:
%
% * asImage - name of the image file to perform algorithm on
% * arSigma - sigma value used for initial Gaussian filtering
% * anWindow - determines size of neighborhood for autocorrelation matrix
% * arThresh - threshold to determine if corner detected
function CornerDetection(asImage, arSigma, anWindow, arThresh)
% Import Image
addpath /Users/jeborah/Desktop/cis/691/hw3/code

lanImage = imread(asImage);

% convert image to greyscale if necessary
if (size(lanImage, 3) ~= 1)
    lanImage = rgb2gray(lanImage);
end
%% Smooth Image and Obtain Gradients
% imgaussfilt performs a convolution with a 5x5 Gaussian kernel
larSmoothedImage = imgaussfilt(lanImage, arSigma);

% Compute the horizontal and vertical gradients
lanGradX = [-1, 0, 1; -1, 0, 1; -1, 0, 1];
lanGradY = lanGradX';
larGradImgX = conv2(larSmoothedImage, lanGradX, 'same');
larGradImgY = conv2(larSmoothedImage, lanGradY, 'same');

% Get the square and product of the components
larGradImgX_2 = larGradImgX.^2;
larGradImgY_2 = larGradImgY.^2;
lanGradImgXY = larGradImgX.*larGradImgY;

% Number of rows and columns
lnNumRows = size(larSmoothedImage,1);
lnNumCols = size(larSmoothedImage,2);

%% Detect All Potential Corners
% Matrix to hold the eigenvalues and coordinates
larL = [];
% Get eigenvalue for each point
for i=1+anWindow:lnNumRows-anWindow
    for j=1+anWindow:lnNumCols-anWindow
        lnSquaredX = ...
            sum(sum(larGradImgX_2(i-anWindow:i+anWindow,...
                j-anWindow:j+anWindow)));
        lnSquaredY = ...
            sum(sum(larGradImgY_2(i-anWindow:i+anWindow,...
                j-anWindow:j+anWindow)));
        lnSquaredProduct= ...
            sum(sum(lanGradImgXY(i-anWindow:i+anWindow,...
                j-anWindow:j+anWindow)));
        
        % Build matrix C 
        larC =[lnSquaredX lnSquaredProduct;...
            lnSquaredProduct lnSquaredY];
        % Store eigenvalues in matrix L
        lnEig = min(eig(larC));
        if lnEig > arThresh
            larL = [larL; [lnEig i j]];
        end
    end
end

%% Filter Weaker Adjacent Corner Detections
% sort the eigenvalues
larL = sortrows(larL, 'descend');

% Keep only the largest eigenvalues in each neighborhood
lnA = 1;
while (lnA ~= length(larL))
    lnB = lnA + 1;
    while (lnB ~= length(larL))
        if (abs(larL(lnA,2) - larL(lnB,2)) <= anWindow && ...
                abs(larL(lnA,3) - larL(lnB,3)) <= anWindow)
            larL(lnB,:) = [];
        else
            lnB = lnB + 1;
        end
    end
    lnA = lnA + 1;
end

%% Plot Corners on Image
imshow(lanImage);
hold on 
for i = 1:length(larL)
    plot(larL(i,3),larL(i,2),'r*', 'MarkerSize', 1);
end
end

