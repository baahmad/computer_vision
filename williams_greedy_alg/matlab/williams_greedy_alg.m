%% Williams/Shah Greedy Algorithm
% Bilal Ahmad 3/8/20

addpath /Users/jeborah/Desktop/

% Import Image
addpath /Users/jeborah/Desktop/cis/691/hw4/code
addpath /Users/jeborah/Desktop/cis/691/hw4/code/Images1through8
addpath /Users/jeborah/Desktop/cis/691/hw4/code/Sequence1
addpath /Users/jeborah/Desktop/cis/691/hw4/code/Sequence2

lanImage = imread("Image2.jpg");
arSigma = 0.5;
lnSize = 1;
lnAlpha = 0;
anBeta = 0;
lrGamma = 1.2;
lrThresh1 = .25;
lrThresh2 = 100;

% convert image to greyscale if necessary
if (size(lanImage, 3) ~= 1)
    lanImage = rgb2gray(lanImage);
end
%% Smooth Image and Obtain Gradient Magnitudes
% imgaussfilt performs a convolution with a 5x5 Gaussian kernel
larSmoothedImage = imgaussfilt(lanImage, arSigma);
larAdjImage = imgradient(larSmoothedImage);

lnNumRows = size(larSmoothedImage,1);
lnNumCols = size(larSmoothedImage,2);

% This matrix will hold the gradient magnitude for each point
larMagGradients = zeros(lnNumRows, lnNumCols);

% Calculate the gradient magnitude for each point
for i = 1+lnSize:lnNumRows-lnSize
    for j = 1+lnSize:lnNumCols-lnSize        
        lnMaxGrad = max(max(larAdjImage(i-lnSize:i+lnSize,...
            (j-lnSize):(j+lnSize))));
        lnMinGrad = min(min(larAdjImage(i-lnSize:i+lnSize,...
            j-lnSize:j+lnSize)));
        
        if (lnMaxGrad - lnMinGrad < 5)
            lnMinGrad = lnMaxGrad - 5;
        end
        larMagGradients(i,j) = (lnMinGrad-larAdjImage(i,j)) / ...
            (lnMaxGrad - lnMinGrad);
    end
end
clear lnMaxGrad lnMinGrad

%% Obtain initial contour points from user
% prompt user to select contours
fprintf("Select initial contour points\n");
imshow(larSmoothedImage)
[lanSelectX,lanSelectY] = getline;

% add extra points if necessary
lanPoints = [];
lnNumSel = length(lanSelectX);

% % use the line equation between points to figure out where
% % the interpolated added ponts go
for i =1:lnNumSel
    if (i == lnNumSel)
        lnNextX = lanSelectX(1);
        lnNextY = lanSelectY(1);
    else
        lnNextX = lanSelectX(i+1);
        lnNextY = lanSelectY(i+1);
    end
    
    lanPoints = [lanPoints ; lanSelectX(i) lanSelectY(i)];
    
    % distance between points
    lrDist = sqrt((lanSelectX(i)-lnNextX)^2 + ...
        (lanSelectY(i)-lnNextY)^2);
    
    % number of points that need to be added
    lnNumAdd = max(0, round(lrDist/5)-1);
    
    if (lnNumAdd > 0)
        for k=1:lnNumAdd
            lrDisRat = (5*k)/lrDist;
            lanPoints = [lanPoints ; ...
                (1-lrDisRat)*lanSelectX(i) + lnNextX*lrDisRat ...
                (1-lrDisRat)*lanSelectY(i) + lnNextY*lrDisRat]; 
        end
    end
end
clear lanSelectX lanSelectY lrDisRat lrDist

% update the count of selected elements
lanPoints = round(lanPoints);
%lanPoints = lanDogma;
lnNumSel = length(lanPoints);
larBetas = ones(1,lnNumSel) * anBeta;
%% Implement the greedy algorithm by Williams and Shah
lnPointsMoved = realmax;
lnCount = 0;
lbCorner = false;
while lnPointsMoved > 4
    % Start by finding average distance between points
    lnPointsMoved = 0;
    lrDistSum = 0.0;
    for j=1:lnNumSel
        if j == lnNumSel
            lrDistSum = lrDistSum + ...
                sqrt((lanPoints(j,1)-lanPoints(1,1))^2 + ...
                ((lanPoints(j,2)-lanPoints(1,2))^2));
        else
            lrDistSum = lrDistSum + ...
                sqrt((lanPoints(j,1)-lanPoints(j+1,1))^2 + ...
                ((lanPoints(j,2)-lanPoints(j+1,2))^2));
        end
    end
    lrAvgDist = lrDistSum/lnNumSel;
    
    
    for i = 1:lnNumSel+1
        % adjust the initial point again after all other points
        % have been adjusted
        if i > lnNumSel
            i = 1;
        end
        % adjust prev/next points for when loop is on first/last points
        if i == 1
            lanPrevPoint = lanPoints(lnNumSel,:);
        else
            lanPrevPoint = lanPoints(i-1,:);
        end
        
        if i == lnNumSel
            lanNextPoint = lanPoints(1,:);
        else
            lanNextPoint = lanPoints(i+1,:);
        end
        
        % Obtain continuity terms 
        larContTerm = ones(2*lnSize+1, 2*lnSize+1);
        b = 1;
        for j=-lnSize:lnSize
            a = 1;
            for k = -lnSize:lnSize
                larContTerm(a,b) = abs(lrAvgDist - ...
                    sqrt((lanPoints(i,1)+j-lanPrevPoint(1,1))^2 + ...
                    ((lanPoints(i,2)+k-lanPrevPoint(1,2))^2)));
                a = a + 1;
            end
            b = b + 1;
        end
        larContTerm = larContTerm./max(max(larContTerm));
        
        % Obtain curvature terms
        larCurvTerm = ones(2*lnSize+1 ,2*lnSize+1);
        b = 1;
        for j=-lnSize:lnSize
            a = 1;
            for k = -lnSize:lnSize
                larCurvTerm(a,b) = (lanPrevPoint(1,1) - ...
                        2*(lanPoints(i,1)+j) + ...
                        lanNextPoint(1,1))^2 + ...
                    (lanPrevPoint(1,2) - ...
                        2*(lanPoints(i,2)+k) + ...
                        lanNextPoint(1,2))^2;
                a = a + 1;
            end
            b = b + 1;
        end
        larCurvTerm = larCurvTerm./max(max(larCurvTerm));
        
        % Obtain image energy term
        larImageTerm = ones(2*lnSize+1 ,2*lnSize+1);
        b = 1;
        for j=-lnSize:lnSize
            a = 1;
            for k = -lnSize:lnSize
                lnRow = lanPoints(i,2);
                lnCol = lanPoints(i,1);
                
                larImageTerm(a,b) = larAdjImage(...
                    lanPoints(i,2)+k, lanPoints(i,1)+j);
                a = a + 1;
            end
            b = b + 1;
        end
        
        
        lnMaxGrad = max(max(larImageTerm(1:(2*lnSize+1),...
            1:(2*lnSize+1))));
        lnMinGrad = min(min(larImageTerm(1:(2*lnSize+1),...
            1:(2*lnSize+1))));
        
        if (lnMaxGrad - lnMinGrad < 5)
            lnMinGrad = lnMaxGrad - 5;
        end
        
        for a=1:(2*lnSize+1)
            for b = 1:(2*lnSize+1)
                larImageTerm(a,b) = (lnMinGrad-larImageTerm(a,b))/...
                    (lnMaxGrad-lnMinGrad);
            end
        end
        
        % Determine where to move point
        lanNewPoint = [];
        lrBeta = larBetas(i);
        lrMinEnergy = lnAlpha * larContTerm(lnSize+1,lnSize+1) + ...
                    lrBeta * larCurvTerm(lnSize+1,lnSize+1) + ...
                    lrGamma * larImageTerm(lnSize+1,lnSize+1);
        lrEngTot = realmax;
        b = 1;
        for j=-lnSize:lnSize
            a = 1;
            for k = -lnSize:lnSize
                lrEngTot = lnAlpha * larContTerm(a,b) + ...
                    lrBeta * larCurvTerm(a,b) + ...
                    lrGamma * larImageTerm(a,b);
                if lrEngTot < lrMinEnergy
                    lrMinEnergy = lrEngTot;
                    lanNewPoint = [lanPoints(i,1)+j lanPoints(i,2)+k];
                end
                a = a + 1;
            end
            b = b + 1;
        end
        if ~isempty(lanNewPoint) && ...
                (lanPoints(i,1) ~= lanNewPoint(1,1) || ...
                lanPoints(i,2) ~= lanNewPoint(1,2))
            lnPointsMoved = lnPointsMoved + 1;
            lanPoints(i,1) = lanNewPoint(1,1);
            lanPoints(i,2) = lanNewPoint(1,2);
            clear lanNewPoint;
        end
    end
    
    % Allow for corners
    % adjust prev/next points for when loop is on first/last points
    larCs = zeros(1,lnNumSel);
    for i=1:lnNumSel
        if i == 1
            lanPrevPoint = lanPoints(lnNumSel,:);
        else
            lanPrevPoint = lanPoints(i-1,:);
        end
        
        if i == lnNumSel
            lanNextPoint = lanPoints(1,:);
        else
            lanNextPoint = lanPoints(i+1,:);
        end
        lanU = [(lanPoints(i,1) - lanPrevPoint(1,1)) ...
            (lanPoints(i,2)-lanPrevPoint(1,2))];
        lrMagU = sqrt(lanU(1,1)^2 + lanU(1,2)^2);
        lanU = lanU/lrMagU;
        lanNextU = [(lanNextPoint(1,1) - lanPoints(i,1)) ...
            (lanNextPoint(1,2) - lanPoints(i,2))];
        lrNextMagU = sqrt(lanNextU(1,1)^2 + lanNextU(1,2)^2);
        lanNextU = lanNextU/lrNextMagU;
        
        larCs(i) = (lanU(1,1) - lanNextU(1,1))^2 + ...
            (lanU(1,2) - lanNextU(1,2))^2;     
    end
    
    for i=1:lnNumSel
        if i == 1
            lrPrevC = larCs(lnNumSel);
        else
            lrPrevC = larCs(i-1);
        end
        
        if i == lnNumSel
            lrNextC = larCs(1);
        else
            lrNextC = larCs(i+1);
        end
        lrMagV = sqrt(lanPoints(i,1)^2 + lanPoints(1,2)^2);
        if larCs(i) > lrPrevC && larCs(i) > lrNextC && ...
                larCs(i) > lrThresh1 && ...
                lrMagV > lrThresh2
            larBetas(i) = 0;
        else
            larBetas(i) = anBeta;
        end
    end
end

