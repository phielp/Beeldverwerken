function beeldverwerken
    pixelValue('kids.tif',127,50,'nearest')
end
function color = pixelValue(image,x,y,method)
    im = imread(image);
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