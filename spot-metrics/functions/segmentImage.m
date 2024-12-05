function [nucleus, spots, bg] = segmentImage(img)
%SEGMENTIMAGE Create nucleus and spot masks using Otsu's method

    % perform bilateral filtering
    filtimg = bfilter2(im2double(img),5,30);

    t = multithresh(filtimg, 2);      % Determine threshold using Otsu's method
    
    % Segment nucleus based on thresholding
    class = imquantize(filtimg, t);
    
    % find spots
    spots = class;
    spots(spots ~= 3) = 0;
    
    otsu = class;
    otsu(otsu ~= 2) = 0;
    otsu = logical(otsu);
    otsu = imfill(otsu,'holes');
    bg = ~otsu;
    % select candidate based on position
    center = false(size(img));
    center(round(size(img,1)/2),round(size(img,2)/2)) = 1;
    center = imdilate(center,strel('disk',50));
    [centercolumn, centerrow] = find(center);
    otsu = bwselect(otsu,centercolumn,centerrow,4);
    
    
    % Segment nucleus based on edge-detection
    
    % automatically detect edge thresholds
    [~,t] = edge(img,'canny');
    
    % perform canny edge detection with bias for low frequency edges
    canny = edge(img,'canny',1.5*t);
    
    % dilate image to close gaps
    canny = imdilate(canny,strel('disk',15));
    
    % fill holes in the detection
    canny = imfill(canny, 'holes');
    
    % remove areas too small to be a nucleus
    canny = bwareaopen(canny, nnz(spots));
    
    % select candidate based on position
    canny = bwselect(canny,centercolumn,centerrow,4);
    
    % take the segmentation to be the areas of agreement
    nucleus = canny & otsu;
    
    % select final segmentation based on position
    nucleus = bwselect(nucleus, centercolumn, centerrow, 4);
    
    % prune segmentation
    center = round(size(nucleus)/2);
    nucleus = bwprune(nucleus, center(1),center(2));
    
    % cut spot detections to within relevant nucleus
    spots = spots & nucleus;
    
    % Throw error if segmentation failed
    if nnz(nucleus) == 0
        error('Segmentation fault');
    end

end

