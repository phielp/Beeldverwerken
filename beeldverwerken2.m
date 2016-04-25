function beeldverwerken2
    %% Theory question part 1
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
    % Expecting: 1
    sum(sum(Gauss(10))) % ans = 1.0000
    sum(sum(Gauss(100)))% ans = 1.0000
    
    mesh(Gauss(3))
end

function G = Gauss(S)
    % create appropriate ranges for x and y
    sigma = S;
    M = 2*sigma;  
    N = 2*sigma;
    x = -M : M;
    y = -N : N;
    % create a sampling grid
    [X, Y] = meshgrid(x,y);
    G=exp(-X.^2/(2*sigma^2)-Y.^2/(2*sigma^2));
    G=G./sum(G(:));
    %G = 1/((sigma*sqrt(2*pi)^2))*exp(-(X.^2+Y.^2)/(2*sigma^2));
end
