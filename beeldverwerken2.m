% LAB 2, 25-4-2016
%
% Authors:
%   Steven de Weille, 10606750
%   Philip Bouman, 10668667

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

    %% Question 3
    % 3.2 Derivatives
    
    % given code
    x = -100:100;
    y = -100:100;
    [X, Y] = meshgrid(x, y);
    A = 1; B = 2; V = 6*pi/201; W = 4*pi/201;
    
    % x and y derivative 
    Fx = A * cos(V * X);
    Fy = -W * B * sin(W * Y);
    
    % show
    figure;
    subplot(1,2,1);
    imshow(Fx, [], 'xData', x, 'yData', y);
    subplot(1,2,2);
    imshow(Fy, [], 'xData', x, 'yData', y);
    
    clear;
    
    % 3.3 Gradient vectors
    
    % given code
    x = -100:100;
    y = -100:100;
    [X, Y] = meshgrid(x, y);
    A = 1; B = 2; V = 6*pi/201; W = 4*pi/201;
    F = A * sin(V*X) + B * cos(W*Y);
    
    % plot arrows
    xx = -100:10:100;
    yy = -100:10:100;
    [XX, YY] = meshgrid(xx, yy);

    % x and y derivative (see 3.1)
    Fx = A*V*cos(V*XX);
    Fy = -B*W*sin(W*YY);

    % show
    figure;
    imshow(F, [], 'xData', x, 'yData', y);
    hold on;
    quiver(XX, YY, Fx, Fy, 'r');
    hold off;
    
    % 3.4 Gradient vectors (follow-up)
    
    % given code
    x = -100:100;
    y = -100:100;
    [X, Y] = meshgrid(x, y);
    A = 1; B = 2; V = 6*pi/201; W = 4*pi/201;
    F = A * sin(V*X) + B * cos(W*Y);
    
    % plot arrows
    xx = -100:10:100;
    yy = -100:10:100;
    [XX, YY] = meshgrid(xx, yy);
    
    % get gradient vectors
    Fx = gD(F, 1, 1, 0);
    Fy = gD(F, 1, 0, 1);
    Fx1 = Fx(xx + 101, yy + 101);
    Fy1 = Fy(xx + 101, yy + 101);
    
    % show
    figure;
    imshow(F, [], 'xData', x, 'yData', y);
    hold on;
    quiver(XX, YY, Fx1, Fy1, 'r');
    hold off;
    
    % 3.5 Gradient vectors on rotated image
    
    % given code
    x = -100:100;
    y = -100:100;
    [X, Y] = meshgrid(x, y);
    A = 1; B = 2; V = 6*pi/201; W = 4*pi/201;
    F = A * sin(V*X) + B * cos(W*Y);
    
    % plot arrows
    xx = -100:10:100;
    yy = -100:10:100;
    [XX, YY] = meshgrid(xx, yy);
    
    % rotate
    R = rotateImage(F, 6/pi, 'linear');
    
    % get gradient vectors
    Fx = gD(F, 1, 1, 0);
    Fy = gD(F, 1, 0, 1);
    Fx = Fx(xx + 101, yy + 101);
    Fy = Fy(xx + 101, yy + 101);
    
    % show
    figure;
    subplot(1, 2, 1);
    title('Rotated with gradient vectors');
    imshow(F, [], 'xData', x, 'yData', y);
    hold on;
    quiver(XX, YY, Fx, Fy, 'r');
    hold off;
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

% 2.8
function [ G ] = gauss1(sigma)

    M = abs(ceil(2.5 * sigma));
    
    sd = 2 * sigma^2;
    
    x = linspace(ceil(-M/2), floor(M/2), M);
    G = exp(-x.^2/sd); 
    
    G = G ./ sum(G(:));
end    

% 2.9
function [ F ] = gD(F, sigma, xorder, yorder)

    % calculates gaussian derivative
    G_x = gauss1(sigma);
    
    M = abs(ceil(2.5 * sigma));
    x = linspace(ceil(-M/2), floor(M/2), M);
    dx = -(x./sigma^2) .* G_x;
    
    if (xorder == 1 && yorder == 0)
        % derivative to x
        F = imfilter(F, dx, 'conv', 'replicate');
    elseif (yorder == 1 && xorder == 0)
        % derivative to y
        F = imfilter(F, dx', 'conv', 'replicate');
    elseif (xorder == 2 && yorder == 0)
        % 2 times to x
        F = imfilter(F, dx, 'conv', 'replicate');
        F = imfilter(F, dx, 'conv', 'replicate');
    elseif (yorder == 2 && xorder == 0)
        % 2 times to y
        F = imfilter(F, dx', 'conv', 'replicate');
        F = imfilter(F, dx', 'conv', 'replicate');
    elseif (yorder == 1 && xorder == 1)
        % derivate to x
        F = imfilter(F, dx, 'conv', 'replicate');
        % derivative to y
        F = imfilter(F, dx', 'conv', 'replicate');
    elseif (xorder == 2 && yorder == 1)
        % 2 times to x
        F = imfilter(F, dx, 'conv', 'replicate');
        F = imfilter(F, dx, 'conv', 'replicate');
        % 1 time to y
        F = imfilter(F, dx', 'conv', 'replicate');
    elseif (yorder == 2 && xorder == 1)
        % 2 times to y
        F = imfilter(F, dx', 'conv', 'replicate');
        F = imfilter(F, dx', 'conv', 'replicate');
        % 1 time to x
        F = imfilter(F, dx, 'conv', 'replicate');
    elseif (yorder == 2 && xorder == 2)
        % 2 times to x
        F = imfilter(F, dx, 'conv', 'replicate');
        F = imfilter(F, dx, 'conv', 'replicate');
        % 2 times to y
        F = imfilter(F, dx', 'conv', 'replicate');
        F = imfilter(F, dx', 'conv', 'replicate');
    else
        assert(0 == 1, 'unsupported');
    end
end   

% Rotation (assignment 1)
function rotatedImage = rotateImage(image, angle, method)
    angle = degtorad(angle);
    [imy, imx] = size(image);
    t1 = imx/2;
    t2 = imy/2;
    % rotation matrix
    R = [cos(angle), -sin(angle), t1; 
         sin(angle), cos(angle), t2; 
         0, 0, 1];
    
    T = [1,0,-t1; 0,1,-t2; 0,0,1];
    [X, Y] = meshgrid(1:imx,1:imy);
    Indices = [X(:) Y(:)]';
    Indices = [Indices; ones(1,length(Indices))];
    Rotated = R*T*Indices;

    for k=1:length(Rotated)
        if Rotated(1,k)<= imx && Rotated(2,k) <= imy
            color(k) = pixelValue(image, Rotated(2,k), Rotated(1,k),method);
        else
            color(k) = 0;
        end
    end
    rotatedImage = reshape(color,imx,imy);
end

% Interpolation (assignment 1)
function color = pixelValue(image,x,y,method)
    im = image;
    [height,width] = size(im);
    if inImage(x,y,width,height)
        switch(method)
            % nearest-neighbour interpolation
            case 'nearest'
                nearestX = floor(x+0.5);
                nearestY = floor(y+0.5);
                if nearestX < 1
                    nearestX = 1;
                end
                if nearestY < 1
                    nearestY = 1;
                end
                if nearestX > width
                    nearestX = width;
                end
                if nearestY > height
                    nearestY = height;
                end
                color = im(nearestX,nearestY);
                return
            % bilinear interpolation    
            case 'linear'
                k = floor(x);
                l = floor(y);              
                a = x-k;
                b = y-l;
                if k < 1
                    k = 1;
                end
                if l < 1
                    l = 1;
                end
                color = (1-a)*(1-b)*im(k,l)...
                +(1-a)*b*im(k,l+1)...
                +a*b*im(k+1,l+1)...
                +a*(1-b)*im(k+1,l);               
        end        
    else 
        color = 1;
    end
end

% image check (assignment 1)
function boolean = inImage(x1,y1,width,height)
    boolean = ((x1 <= width) && (y1 <= height) && (height > 1) && (width > 1));
end
