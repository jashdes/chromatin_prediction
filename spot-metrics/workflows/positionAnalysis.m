function [pos, sigma,p,normp,ratio, intensity] = positionAnalysis(greenimage,redimage,template, scale)
%POSITIONANALYSIS Performs position analysis on a grid image

%   image:      image containing the grid
%   template:   template of the grid to be found
%   scale:      image scale in micrometers per pixel
%   pos:        n x 2 matrix with each row containing one spot's position

%   p:          distance to the border in pixels
%   normp:      p normalized by dividing by squareroot of nucleus area
%   ratio:      c / (c + p) where c is distance to nucleus center in pixels

% crop input images to grid nucleus neighborhood
[greenimage, redimage] = cropFullFrame(greenimage, {redimage});

% extract redimage from cell array
redimage = redimage{1};

% find the template grid within the image
[grid,xoffset,yoffset] = findGrid(greenimage,template);

% find spot positions within the grid
[gridpos, sigma] = findSpots(grid,template);

% calculate image positions from grid-relative positions
pos = gridpos + [yoffset,xoffset];

% segment the image into nucleus and background
[nucleus, ~, background] = segmentImage(greenimage);

% perform intensity analysis
intensity = getSpotIntensities(redimage,pos,background,nucleus);

% find center of the nucleus
[xcenter,ycenter] = findInnermost(nucleus,0.1);

% perform position analysis
[p, normp, ratio] = getPositionMetrics(pos,[xcenter,ycenter],nucleus);

% scale spatial results
pos = pos * scale;
sigma = sigma * scale;
p = p * scale;

    
end

