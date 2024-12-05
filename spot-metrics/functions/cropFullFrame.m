function [cropped, samecrop] = cropFullFrame(image, samecrop, cropsize)
%CROPFULLFRAME Crops full frame images to ROI containing grid nucleus

    %   image:      image containing the nucleus with grid pattern
    %   samecrop:   image stack cell array to receive the same crop as image
    %   cropsize:   side-length of the square of resulting crop
    
    if(nargin < 3)
        cropsize = 550;
    end
    
    % bias nucleus search towards the center of the image
    bias = zeros(size(image));
    bias(size(image,1)/2,size(image,2)/2) = 1;
    bias = bwdist(bias);
    bias = bias ./ max(max(bias));
    bias = 1 - bias;
    bias = bias .*bias;
    biased = double(image) .* bias;
    
    % find grid nucleus by assuming it will be brightest
    [y,x] = find(biased == max(max(biased)));
    
    % crop area around it
    cropoffset = [x-cropsize/2,y-cropsize/2,cropsize,cropsize];
    cropped = imcrop(image,cropoffset);
    
    % crop stack if it exists
    if(nargin >= 2)
        for k = 1:size(samecrop,1)
            samecrop{k} = imcrop(samecrop{k},cropoffset);
        end
    end
    
   
end
