function [res] = batchCrop(path, args)
%BATCHCROP Batch crop green/red images with grids in their center

    % path:             path to an image stack to be cropped
    % args.cropsize:    target length of the cropped square images
    % args.outpath:     path to write result to
    % args.template:    grid template to look for when cropping
    
    
    
    % load image stack
    imstack = loadImageStack(path);
    
    % crop stack based on first image
    if isfield(args,'cropsize')
        % use custom crop size if provided
        [~,croppedstack] = cropFullFrame(imstack{1},imstack, args.cropsize);
    else
        % userwise fallback to default
        [~,croppedstack] = cropFullFrame(imstack{1},imstack);
    end
    
    % fine crop based on grid
    [~,~,~,croppedstack] = findGrid(croppedstack{1},args.template,croppedstack);
    
    % construct output file name
    [~,name,ext] = fileparts(path);
    outpath = fullfile(args.outpath,strcat(name,ext));
    
    % write cropped stack
    saveImageStack(croppedstack, outpath);
    
    % set result to reflect success
    res = "cropped";
    
end
