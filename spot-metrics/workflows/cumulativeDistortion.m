function [len,ang] = cumulativeDistortion(img, template, tpos)
%CUMULATIVEDISTORTION get cumulative distortion accross all frames of img

%   img:        image stack to be analyzed
%   template:   template to use to find spots within image
%   tpos:       optional - tracking output to use instead of image
%   len:        cumulative change in lengths across all frames
%   ang:        cumulative change in angles across all frames

% extract position from image stack
if(nargin < 3)
    tpos = spotTracking(img,template);
end

% Reshape tracking output for distortion analysis
    for k = 1:size(tpos,1)
        pos(k,:,1) = tpos{k}.x;
        pos(k,:,2) = tpos{k}.y;
    end

% Calculate distoritions for each frame
[len,ang] = gridDistortion(pos);

% Calculate cumulative distortions
len = cumsum(len);
ang = cumsum(ang);
    
end