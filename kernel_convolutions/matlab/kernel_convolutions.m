%% Kernel Convoltions
% Bilal Ahmad

%% Import Images
addpath /Users/jeborah/Desktop/

lanNoisyImage1 = imread('NoisyImage1.jpg');
lanNoisyImage2 = imread('NoisyImage2.jpg');

%% Apply Mean Filter

% Make the kernels
lanKernel1 = ones(3,3)/9;
lanKernel2 = ones(9,9)/81;

% Convolve keeping size of the images
larK1_Image1 = conv2(lanNoisyImage1, lanKernel1, 'same'); 
larK2_Image1 = conv2(lanNoisyImage1, lanKernel2, 'same'); 
larK1_Image2 = conv2(lanNoisyImage2, lanKernel1, 'same'); 
larK2_Image2 = conv2(lanNoisyImage2, lanKernel2, 'same'); 

% Save the comparison plots
subplot (1,3,1)
imshow(lanNoisyImage1);
title('Original Image');
subplot(1,3,2)
imshow(larK1_Image1, []);
title('3x3 Mean Kernel');
subplot(1,3,3)
imshow(larK2_Image1, []);
title('9x9 Mean Kernel');
saveas(gcf, 'Mean_I1', 'eps')

subplot (1,3,1)
imshow(lanNoisyImage2);
title('Original Image');
subplot(1,3,2)
imshow(larK1_Image2, []);
title('3x3 Mean Kernel');
subplot(1,3,3)
imshow(larK2_Image2, []);
title('9x9 Mean Kernel');
saveas(gcf, 'Mean_I2', 'eps')

%% Apply Gaussian Filter

% Make two 9x9 kernels with different sigma values
lnSigma1 = 1;
lnSigma2 = 3;
larKernel1 = ones(9,9);
larKernel2 = ones(9,9);

lrSum = 0.0;
for i = 1:9
    for j = 1:9
        % for 9x9 kernel the center is (5,5)
        sq_dist = (i-5)^2 + (j-5)^2; 
        larKernel1(i,j) = exp(-1 * (sq_dist)/ (2 * lnSigma1^2));
        lrSum = lrSum + larKernel1(i,j);
    end
end
larKernel1 = larKernel1/lrSum;

lrSum = 0.0;
for i = 1:9
    for j = 1:9
        % for 9x9 kernel the center is (5,5)
        sq_dist = (i-5)^2 + (j-5)^2; 
        larKernel2(i,j) = exp(-1 * (sq_dist)/ (2 * lnSigma2^2));
        lrSum = lrSum + larKernel2(i,j);
    end
end
larKernel2 = larKernel2/lrSum;

% Convolve keeping size of the images
larK1_Image1 = conv2(lanNoisyImage1, larKernel1, 'same'); 
larK2_Image1 = conv2(lanNoisyImage1, larKernel2, 'same'); 
larK1_Image2 = conv2(lanNoisyImage2, larKernel1, 'same'); 
larK2_Image2 = conv2(lanNoisyImage2, larKernel2, 'same'); 

% Save the comparison plots
subplot (1,3,1)
imshow(lanNoisyImage1);
title('Original Image');
subplot(1,3,2)
imshow(larK1_Image1, []);
title('9x9 Kernel, Sigma = 1');
subplot(1,3,3)
imshow(larK2_Image1, []);
title('9x9 Kernel, Sigma = 3');
saveas(gcf, 'Gaussian_I1', 'eps')

subplot (1,3,1)
imshow(lanNoisyImage2);
title('Original Image');
subplot(1,3,2)
imshow(larK1_Image2, []);
title('9x9 Kernel, Sigma = 1');
subplot(1,3,3)
imshow(larK2_Image2, []);
title('9x9 Kernel, Sigma = 3');
saveas(gcf, 'Gaussian_I2', 'eps')

%% Apply Median Filter

% Use medfilt2 to perform convolution
larK1_Image1 = medfilt2(lanNoisyImage1, [3 3]); 
larK2_Image1 = medfilt2(lanNoisyImage1, [9 9]); 
larK1_Image2 = medfilt2(lanNoisyImage2, [3 3]); 
larK2_Image2 = medfilt2(lanNoisyImage2, [9 9]); 

% Save the comparison plots
subplot (1,3,1)
imshow(lanNoisyImage1);
title('Original Image');
subplot(1,3,2)
imshow(larK1_Image1, []);
title('3x3 Median Kernel');
subplot(1,3,3)
imshow(larK2_Image1, []);
title('9x9 Median Kernel');
saveas(gcf, 'Median_I1', 'eps')

subplot (1,3,1)
imshow(lanNoisyImage2);
title('Original Image');
subplot(1,3,2)
imshow(larK1_Image2, []);
title('3x3 Median Kernel');
subplot(1,3,3)
imshow(larK2_Image2, []);
title('9x9 Median Kernel');
saveas(gcf, 'Median_I2', 'eps')