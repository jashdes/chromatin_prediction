function [] = showVideoPair(a,b,fps,outpath)
%SHOWVIDEOPAIR ShowVideo for two image sequences to be compared

    % a:    first image sequence
    % b:    second image sequence
    % fps:  frames per second
    % outpath:  path to which video will be written
    
    % prepare stack to be written
    diff = cell(size(a,1));
    
    % check sequence legnths
    if size(a,1) == size(b,1)
        
        % prepare overlaid frames
        for f = 1 : size(a,1)
            aframe = imadjust(a{f});
            bframe = imadjust(b{f});
            diff{f} = imfuse(aframe,bframe);
        end
        
        showVideo(diff,fps,outpath);
            
    else
        error("Image sequences have different length");
    end
    
end
