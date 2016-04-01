function beeldverwerken
    pixelValue('kids.tif',253,50,'nearest')
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
        end        
    else 
        color = 0;
    end
end

function boolean = inImage(x1,y1,width,height)
    boolean = ((x1 < width) && (y1 < height) && (height > 1) && (width > 1));
end
