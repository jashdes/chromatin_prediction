function [] = timeLapse(path, outpath, fps, motion, register)
%TIMELAPSE Create a timelapse from an image stack

    % path:     path to .tiff image stack representing sequence
    % outpath:  destination to which timelapse will be written
    % fps:      frames per second of resulting timelapse
    % motion:   if set, shows diff to first frame for every frame
    % register: if set, stack will be registered
    
    % load image stack from path
    imstack = loadImageStack(path);
    
    % crop stack based on first image
    [~,imstack] = cropFullFrame(imstack{1},imstack);
    
    % register
    if nargin >= 5 && register
        imstack = stackReg(imstack);
    end
        
    % write video
    showVideo(imstack,fps,strcat(outpath,'.avi'),motion);

end

