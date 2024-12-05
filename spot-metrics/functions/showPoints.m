function [] = showPoints(img, x, y)
%SHOWPOINTS Displays points on top of an image
    imshow(img);
    hold on;
    plot(y, x, 'kd', 'LineWidth', 2, 'MarkerSize', 10);
end

