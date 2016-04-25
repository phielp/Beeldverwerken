
function beeldverwerken
%% Question 1.20
% Read image
im = imread('cameraman.jpg');
im = im2double(im);
im = rgb2gray(im);

% Motionblur (shift 5 px to the right)
mblur = 1/6 * [1 1 1 1 1 1];
% Intensity increase by 3
intens = 3;
% Average 3x3 neighbourhood
average = 1/9 *[1 1 1; 1 1 1; 1 1 1];

% Show results
figure;
subplot(2, 2, 1);
imshow(imfilter(im, mblur, 'conv', 'replicate'));
title('Motion blur with 5px to right');
subplot(2, 2, 2);
imshow(imfilter(im, intens, 'conv', 'replicate'));
title('Increased intensity (+3)');
subplot(2, 2, 3);
imshow(imfilter(im, average, 'conv', 'replicate'));
title('Average from 3x3 neighbourhood');

%% Implementation of Gaussian derivatives
% create appropriate ranges for x and y
x = -M : M ;
y = -N : N ;
% create a sampling grid
[X , Y ] = meshgrid (x , y );
% determine the scale
sigma = S ;
% calculate the Gaussian function
G = gauss(sigma);

function G = gauss(sigma)
    G = 0;
end
