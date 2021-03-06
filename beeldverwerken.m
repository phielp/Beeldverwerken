% LAB 1, 11-4-2016
%
% Authors:
%   Steven de Weille, 10606750
%   Philip Bouman, 10668667

%% Main
% Questions 2 - 7
% run program command: beeldverwerken
% all plots and figures will show
function beeldverwerken
%% Question 2
% Two interpolation methods: nearesst neighbour en linear
    im = imread('autumn.tif');
    im = im2double(rgb2gray(im));
    n = 100;
    reset(gcf);
    hold on;
    figure(1);
    plot(profile(im, 100, 100, 120, 120, n, 'linear'), 'b');
    plot(profile(im, 100, 100, 120, 120, n, 'nearest'), 'r');
    hold off;
    
    
%% Question 3
% Rotates an image.
    im = imread('cameraman.tif');
    im = im2double(im);
    degrees = 45;
    rotated = rotateImage(im, degrees, 'nearest');
    figure(2);
    imshow(rotated);


%% Question 4 
% Transforms a predefined parallelogram-shaped cut out
% into a straightend square shape.

    im = imread('cameraman.tif');
    im = im2double(im);
    
    % transformed points
    imSize = size(im);
    v1 = [imSize(1)/2-10,10];
    v2 = [imSize(1)/2-120, 100];
    v3 = [imSize(1)/2-20, 180];
    v4 = v1 + (v3 - v2);
    
    % transform image
    x1 = v1(1);
    y1 = v1(2);
    x2 = v2(1);
    y2 = v2(2);
    x3 = v4(1);
    y3 = v4(2);
    transformed = myAffine(im, x1, y1, x2, y2, x3, y3, 256, 256, 'linear');
   
    figure(3);
    subplot(1, 3, 1);
    imshow(transformed);
    subplot(1, 3, 2);
    imshow(im);
    plotParallelogram(x1, y1, x2, y2, v3(1), v3(2));
   
%% Question 5
% Prompts the user to select four points to project.
% For correct results, please click in the following order:
%   -top left corner
%   -bottom left corner
%   -top right corner
%   -bottom right corner

    figure(4);
    im = imread('flyers.png');
    im = imresize(im, 0.6);
    im = im2double(im);
    imshow(im);
    
    % store user selection
    points = ginput(4);
    x1 = points(1, 1);
    y1 = points(1, 2);
    x2 = points(2, 1);
    y2 = points(2, 2);
    x3 = points(3, 1);
    y3 = points(3, 2);
    x4 = points(4, 1);
    y4 = points(4, 2);
    
    projection = myProjection(im, x1, y1, x2, y2, x3, y3, x4, y4, 150, 300, 'linear');
    figure(5);
    imshow(projection);
    
%% Question 7 
% Creates projection matrix

    load calibrationpoints;
    estimatedprojection = estimateProjectionmatrix(xy, XYZ);
    
end

%% FUNCTIONS

% 2) Interpolation
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

% 2) Check if point is inside image 
function boolean = inImage(x1,y1,width,height)
    boolean = ((x1 <= width) && (y1 <= height) && (height > 1) && (width > 1));
end

% 2) Image profile
function line = profile(image, x0, y0, x1, y1, n, method)
    x = linspace(x0, x1, n); 
    y = linspace(y0, y1, n);
    for i = 1:length(x)
        line(i) = pixelValue(image, x(i), y(i), method);
    end
end

% 3) Rotation
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

% 4) Affine Transformation
function r = myAffine(image, x1, y1, x2, y2, x3, y3, M, N, method)
    r = zeros(N, M);
    % create transformation matrix
    A = [0,0,1; N,0,1; 0,M,1]';
    B = [x1, y1; x2, y2; x3, y3]';
    X = B / A;
    
    % transform coordinates
    for xa = 1:M
        for ya = 1:N
            p = X * [ya,xa,1]';
            x = p(1);
            y = p(2);
            r(xa, ya) = pixelValue(image, x, y, method);
        end
    end
end

% 4) plot parallelogram on the current axis
function plotParallelogram(x1, y1, x2, y2, x3, y3)
    hold on;
    plot([x1, x2, x3, x3-x2+x1, x1], [y1, y2, y3, y1-y2+y3, y1], 'y', 'LineWidth', 2);
    text(x1, y1, '1', 'Color', [0, 1, 0], 'FontSize', 18);
    text(x2, y2, '2', 'Color', [0, 1, 0], 'FontSize', 18);
    text(x3, y3, '3', 'Color', [0, 1, 0], 'FontSize', 18);
end

% 5) Create a projection matrix
function projMatrix = createProjectionMatrix(xy, uv)
    % code from assignment
    x = xy(:, 1);
    y = xy(:, 2);
    u = uv(:, 1);
    v = uv(:, 2);
    o = ones(size(x));
    z = zeros(size(x));
    Aoddrows = [x, y, o , z, z, z, -u.*x, -u.*y, -u];
    Aevenrows = [z, z, z, x, y, o, -v.*x, -v.*y, -v];
    A = [Aoddrows; Aevenrows];
    
    % SVD to get projMatrix
    [U,D,V] = svd(A);
    projMatrix = V(:,end);
    
    % reshape into correct dimensions
    projMatrix = reshape(projMatrix, 3, 3);
    projMatrix = projMatrix';
end

% 5) Projection
function projection = myProjection(image, x1, y1, x2, y2, x3, y3, x4, y4, m, n, method)
    % create empty matrix and fill coordinates
    projection = zeros(n, m);
    xy = [[0, 0]; [n, 0]; [0, m]; [n, m]];
    uv = [[x1, y1]; [x2, y2]; [x3,y3]; [x4, y4]];
      
    % loop over projection matrix
    projectionMatrix = createProjectionMatrix(xy, uv);
    for xIndex = 1:m
        for yIndex = 1:n
            
            % projection
            index = [yIndex, xIndex, 1]';
            vec = projectionMatrix * index;
            
            % get real coordinates
            vec = vec / vec(3);
            x = vec(2);
            y = vec(1);
            
            % get pixel values
            projection(yIndex, xIndex) = pixelValue(image, x, y, method);  
        end
    end
end

% 7) estimated projection matrix
function matrix = estimateProjectionmatrix(xy, XYZ)

    x = xy(:, 1);
    y = xy(:, 2);
    X = XYZ(:, 1);
    Y = XYZ(:, 2);
    Z = XYZ(:, 3);
    o = ones(size(x));
    z = zeros(size(x));
    Aoddrows = [X, Y, Z, o, z, z, z, z, -x.*X, -x.*Y, -x.*Z, -x];
    Aevenrows = [z, z, z, z, X, Y, Z, o, -y.*X, -y.*Y, -y.*Z, -y];
    A = [Aoddrows; Aevenrows];
    
    [U,S,V] = svd(A);
    matrix = V(:,end);
    matrix = reshape(matrix, 3, 4);
    matrix = matrix';
end