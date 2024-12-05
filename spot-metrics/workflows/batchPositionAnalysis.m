function [res] = batchPositionAnalysis(path, args)
%BATCHPOSITIONANALYSIS Used with forEachFile to do position analysis on a folder

    % path:             path to an image to be analyzed
    % args.template:    template to be used in analysis
    % args.scale:       micrometers per pixel
    % res:              result struct containing path and position metrics
    
    % load image
    greenimage = imread(path);
    redimage = imread(path,2);
    
    % perform position analysis on the input file
    [pos,sigma,p,normp,ratio, intensity] = ...
        positionAnalysis(greenimage,redimage,args.template, args.scale);
    
    % package results
    
    % position
    res.x = pos(:,1);
    res.y = pos(:,2);
    
    % sigma
    res.sigmaX = sigma(:,1);
    res.sigmaY = sigma(:,2);
    
    % position metrics
    res.p = p;
    res.normp = normp;
    res.ratio = ratio;
    
    % intensity
    res.intensity = intensity;
    
end
