function [stack] = loadImageStack(path)
%LOADIMAGESTACK Load image stack located at path into cell array

    % improved method to load image stack from
    % https://blogs.mathworks.com/steve/2009/04/02/matlab-r2009a-imread-and-multipage-tiffs/
    
    info = imfinfo(path);
    num_images = numel(info);
    stack = cell(num_images,1);
    for f = 1:num_images
        stack{f} = imread(path, f, 'Info', info);
    end

end
