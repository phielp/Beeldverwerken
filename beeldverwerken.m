% Main
function beeldverwerken
    rotateImage('cameraman.tif', 30, 'linear');
end

function interpolationPlot
    a = imread('autumn.tif');
    a = im2double(rgb2gray(a));
    n = 100;
    reset(gcf);
    hold on;
    plot(profile(a, 100, 100, 120, 120, n, 'linear'), 'b');
    plot(profile(a, 100, 100, 120, 120, n, 'nearest'), 'r');
    hold off;
end

% Interpolation
function color = pixelValue(image,x,y,method)
    %im = imread(image);
    im = image;
    [height,width] = size(im);
    if inImage(x,y,width,height)
        switch(method)
            case 'nearest'
                nearestX = floor(x+0.5);
                nearestY = floor(y+0.5);
                color = im(nearestX,nearestY);
                return;
            case 'linear'
                k = floor(x+0.5);
                l = floor(y+0.5);
                a = x-k;
                b = y-l;
                color = (1-a)*(1-b)*im(k,l)...
                +(1-a)*b*im(k,l+1)...
                +a*b*im(k+1,l+1)...
                +a*(1-b)*im(k+1,l);               
        end        
    else 
        color = -1;
    end
end

% Check if point is inside image 
function boolean = inImage(x1,y1,width,height)
    boolean = ((x1 <= width) && (y1 <= height) && (height > 1) && (width > 1));
end

% Profile
function line = profile(image, x0, y0, x1, y1, n, method)
    x = linspace(x0, x1, n); 
    y = linspace(y0, y1, n);
    for i = 1:length(x)
        line(i) = pixelValue(image, x(i), y(i), method);
    end
end

% Rotation
function rotatedImage = rotateImage(image, angle, method)
     
    im = imread(image);
    im = im2double(im);

    % image size
    [image_x, image_y] = size(im);

    % center
    c = [image_x; image_y] / 2;

    % rotation matrix
%     R = [cos(-angle), -sin(-angle), c(1); 
%         sin(-angle), cos(-angle), c(2); 
%         0, 0, 1]

    R = [cos(-angle), -sin(-angle);
        sin(-angle), cos(-angle)];
    
    for x = 1:image_x
        for y = 1:image_y
            
            p = [x; y] - c;
            p = R * p;
            p = p + c;
            
            rotated_image(x, y) = pixelValue(image, p(1), p(2), method);
            
        end
    end
    imshow(rotated_image);
end
    

