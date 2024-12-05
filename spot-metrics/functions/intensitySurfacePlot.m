function [] = intensitySurfacePlot(map)
%INTENSITYSURFACEPLOT surface plot with image intensities as third axis

    [x,y] = meshgrid(1:size(map,2),1:size(map,1));
    figure; 
    mesh(x,y,map);

end

