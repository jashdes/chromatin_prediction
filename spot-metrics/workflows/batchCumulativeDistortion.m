function [res] = batchCumulativeDistortion(path, args)
%BATCHCUMULATIVEDISTORTION CumDist on folders

    % path:             path to an image to be analyzed
    % args.template:    template to be used for spot matching
    % res:              result struct containing final distortions
    
    [~,~,ext] = fileparts(path);
    
    if(strcmp(ext,'.txt'))
        % load positions from text file
        tpos = importTxtPositions(path);
        % perform position analysis on the text file
        [len,ang] = cumulativeDistortion(false, false, tpos);
    else
        % load image
        img = loadImageStack(path);

        % strip stack keeping first and last frame
        img = {img{1};img{end}};

        % perform position analysis on the input file
        [len,ang] = cumulativeDistortion(img,args.template);
    end
    
    % return only the total distortion
    res.len = len(end);
    res.ang = ang(end);
    
end
