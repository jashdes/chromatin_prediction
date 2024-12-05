function [lesioned] = lesionImage(image, circle)
%LESIONIMAGE Lesions input image by averaging out a circular subset

%   image:      image to be lesioned
%   circle:     standard MATLAB circle notation [x,y,radius]
%   lesioned:   lesioned image

    % draw the specified circle into an RGB image
    mask = insertShape(zeros(size(image)),'FilledCircle', circle);
    
    % Strip RGB to grayscale to create a logical mask
    mask = logical(mask(:,:,1));
    
    % Obtain average value of image within the mask
    avg = image;
    avg = mean(mean(avg(mask)));
    
    lesioned = image;
    lesioned(mask) = avg;

end
