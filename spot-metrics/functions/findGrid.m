function [grid,x,y, samecrop] = findGrid(image,template, samecrop)
%registerImages Find the grid in the center of an image

%   image:      image containing the nucleus with the grid in its center
%   template:   template of grid to be found
%   x,y:        grid coordinate offset to offset relative positions

image = im2double(image);
template = im2double(template);

% Default spatial referencing objects
fixedRefObj = imref2d(size(image));
movingRefObj = imref2d(size(template));

% Intensity-based registration
[optimizer, metric] = imregconfig('multimodal');
metric.NumberOfSpatialSamples = 500;
metric.NumberOfHistogramBins = 50;
metric.UseAllPixels = true;
optimizer.GrowthFactor = 1.050000;
optimizer.Epsilon = 1.50000e-06;
optimizer.InitialRadius = 6.25000e-03;
optimizer.MaximumIterations = 500;

% Align centers
fixedCenterXWorld = mean(fixedRefObj.XWorldLimits);
fixedCenterYWorld = mean(fixedRefObj.YWorldLimits);
movingCenterXWorld = mean(movingRefObj.XWorldLimits);
movingCenterYWorld = mean(movingRefObj.YWorldLimits);
translationX = fixedCenterXWorld - movingCenterXWorld;
translationY = fixedCenterYWorld - movingCenterYWorld;

% Coarse alignment
initTform = affine2d();
initTform.T(3,1:2) = [translationX, translationY];

% Apply transformation
tform = imregtform(template,movingRefObj,image,fixedRefObj, ...
    'translation',optimizer,metric,'PyramidLevels',3, ...
    'InitialTransformation',initTform);

% Get local coordinates
[x, y] = transformPointsForward(tform, 0, 0);

% Round coordinates before cropping
x = round(x);
y = round(y);

% Use local coordinates for cropping
cropoffset = [x,y, ...
    size(template,1)-1, size(template,2)-1];
grid = imcrop(image, cropoffset);

% crop stack if given
if nargin >= 3
    for k = 1:size(samecrop,1)
            samecrop{k} = imcrop(samecrop{k},cropoffset);
    end
end

end

