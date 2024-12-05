function [stack] = loadImageStack(path, num_images)
%LOADIMAGESTACK Load image stack located at path into cell array

    % improved method to load image stack from
    % https://blogs.mathworks.com/steve/2009/04/02/matlab-r2009a-imread-and-multipage-tiffs/
    
    
    
    info = imfinfo(path);
    if nargin<2
        num_images = numel(info);
    end
    stack = cell(num_images,1);
    for f = 1:num_images
        stack{f} = imread(path, f, 'Info', info);
    end

end
