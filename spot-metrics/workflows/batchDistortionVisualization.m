function [res,err] = batchDistortionVisualization(path, args)
%BATCHDISTORTIONVISUALIZATION Creates a distortion visualization for images

    % path:     path to image
    % args.out: location to which to write the visualization
    
    % return empty cell arrays
    res = {};
    err = {};
    
    % generate output file name
    [~, name] = fileparts(path);
    p = fullfile(args.out, name);
    
    showDistortion(path,p);
    
end

