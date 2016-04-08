% Main
function beeldverwerken
    %imshow(rotateImage('cameraman.tif', 0, 'linear'));
    im = imread('cameraman.tif');
    im = im2double(im);
    rotated = rotateImage(im, 180, 'nearest');
    imshow(rotated)
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

function color = pixelValue(image,x,y,method)
    %im = imread(image);
    im = image;
    [height,width] = size(im);
    if inImage(x,y,width,height)
        switch(method)
            case 'nearest'
                nearestX = floor(x+0.5);
                nearestY = floor(y+0.5);
                if nearestX < 1
                    nearestX = 1;
                end
                if nearestY < 1
                    nearestY = 1;
                end
                %disp([x,y,nearestX,nearestY])
                if nearestX > width
                    nearestX = width;
                end
                if nearestY > height
                    nearestY = height;
                end
                color = im(nearestX,nearestY);
                return
            case 'linear'
                k = floor(x);
                l = floor(y);              
                a = x-k;
                b = y-l;
                color = (1-a)*(1-b)*im(k,l)...
                +(1-a)*b*im(k,l+1)...
                +a*b*im(k+1,l+1)...
                +a*(1-b)*im(k+1,l);               
        end        
    else 
        %disp('huh')
        color = 1;
    end
end

% Rotation
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
    %disp(Rotated)
    for k=1:length(Rotated)
        if Rotated(1,k)<= imx && Rotated(2,k) <= imy
            color(k) = pixelValue(image, Rotated(2,k), Rotated(1,k),method);
        else
            color(k) = 0;
        end
    end
    rotatedImage = reshape(color,imx,imy);
end



