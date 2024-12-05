function [] = showDistortion(stackpath, outpath)
%SHOWDISTORTION Visualize distortion in an image stack

%   stackpath:  path to image stack 
%   outpath:    path to write resulting image to
%               image will be shown on screen if omitted

        % load image stack information
        info = imfinfo(stackpath);
        num_images = numel(info);
        first = imread(stackpath, 1, 'Info', info);
        last = imread(stackpath, num_images, 'Info', info);
        
        % center images
        [first,last] = cropFullFrame(first,{last});
        last = last{1};
        
        % segment both frames
        firstseg = segmentImage(first);
        lastseg = segmentImage(last);
        
        % generate boundary fused images
        
        % convert to double
        first = im2double(first);
        last = im2double(last);
        
        %  normalize intensity
        first = first / max(first(:));
        last = last / max(last(:));
        
        % fuse registration outlines with normalized images
        firstseg = bwperim(firstseg);
        first = first + firstseg;
        
        lastseg = bwperim(lastseg);
        last = last + lastseg;
        
        % calculate diff
        diff = imfuse(last, first);
        
        if nargin < 2
            imshow(diff);
        else
            imwrite(diff,outpath,'tiff');
        end

end

