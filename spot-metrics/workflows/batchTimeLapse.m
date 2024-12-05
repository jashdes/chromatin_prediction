function [res,err] = batchTimeLapse(path, args)
%BATCHTIMELAPSE Creates a time lapse for image stack at path

    % path:     path to image stack
    % args.fps: frames per second for the timelapse
    % args.outpath: location to which to write the timelapse
    % args.motion:  whether diff to first frame should be shown
    % args.register:    whether stack should be registered
    
    % return empty cell arrays
    res = {};
    err = {};
    
    % generate output file name
    [~, name] = fileparts(path);
    p = fullfile(args.outpath, name);
    
    timeLapse(path,p,args.fps, args.motion, args.register);
    
end

