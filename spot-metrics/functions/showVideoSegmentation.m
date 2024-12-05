function [] = showVideoSegmentation(img,seg, fps, out)
%SHOWVIDEOSEGMENTATION Visualize the result of video segmentation

    %   img:    image stack
    %   seg:    stack of logical matrices containing segmentation
    %   fps:    frame rate of output video
    %   out:    path to which video will be written
    
    % Convert segmentation images to segmentation outlines
    seg = cellfun(@bwperim,seg,'un',0);
    
    % Convert logical images in segmentation stack to double images
    seg = cellfun(@double,seg, 'un',0);
    
    % Lower the max intensity to avoid cropping
    seg = cellfun(@(x) x*0.95,seg,'un',0);
    
    showVideoPair(img,seg,fps,out);
    
    % release memory
    clear(seg);
    
end

