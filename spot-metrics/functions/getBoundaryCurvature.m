function [curvature] = getBoundaryCurvature(img, segmentlength)
%GETBOUNDARYCURVATURE Get curvature an object's boundary from segmentation

    %   img:            logical image representing object and background
    %   segmentlength:  how many pixels are used for fitting a segment
    %   curvature:      curvatures of the fitted segments
    
    %   adapted from https://www.mathworks.com/matlabcentral/answers/164349-how-to-calculate-the-curvature-of-a-boundaries-in-binary-images
    
    % generate boundary from segmentation
    bound = bwboundaries(img,8, 'noholes');
    
    % strip axes
    x = bound{1}(:, 2);
    y = bound{1}(:, 1);
    
    % allocate curvature array
    curvature = NaN(length(bound{1}),1);
    
    % calculate curvatures
    center = floor(segmentlength/2);
    
    % TODO: wrap around beginning and endp
    for k = center+1 : length(x) - center
        
        % get adjacent points
        adjacentX = x(k-center:k+center);
        adjacentY = y(k-center:k+center);
        
        % set selected points within logical image
        adjacent = false(size(img));
        adjacent(adjacentY,adjacentX) = true;
        
        % get centroid and orientation of selected points
        centroid = regionprops(adjacent,'centroid');
        centroid = centroid.Centroid;
        orientation = regionprops(adjacent,'orientation');
        orientation = orientation.Orientation;
        
        % move points to origin
        adjacentX = adjacentX - centroid(1);
        adjacentY = adjacentY - centroid(2);
        
        % rotate points to align with x axis
        R = [cosd(-orientation) -sind(-orientation); ...
            sind(-orientation) cosd(-orientation)];
        rotated = NaN(length(adjacentX),2);
        for r = 1 : length(adjacentX)
            rotated(r, :) = (R * [adjacentX(r); adjacentY(r)]).';
        end
        
        % Get a fit.
        coefficients = polyfit(rotated(:,1), rotated(:,2), 2);
        % Get the curvature
        curvature(k) = abs(coefficients(1)); 
    end
    
end

