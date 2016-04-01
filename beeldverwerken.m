function beeldverwerken
    a = imread('autumn.tif');
    a = im2double(rgb2gray(a));
    plot(profile(a, 100, 100, 120, 120, 100, 'linear'));
end

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

function boolean = inImage(x1,y1,width,height)
    boolean = ((x1 <= width) && (y1 <= height) && (height > 1) && (width > 1));
end

function line = profile(image, x0, y0, x1, y1, n, method)
    % profile of an image along straight line in n points
    x = linspace(x0, x1, n); 
    y = linspace(y0, y1, n);
    for i = 1:length(x)
        line(i) = pixelValue(image, x(i), y(i), method);
    end
end

