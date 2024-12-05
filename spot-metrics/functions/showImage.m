function [] = showImage(path)
%SHOWIMAGE Displays the contrast-adjusted image pointed to by path
    imshow(imadjust(imread(path)));
    [~,n,~] = fileparts(path);
    title(n, 'Interpreter', 'none');
end
