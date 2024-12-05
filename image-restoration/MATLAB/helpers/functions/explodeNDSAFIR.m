function [] = explodeNDSAFIR(stack, path, name)
%EXPLODESTACK writes argument stack frames to folder

    stack = loadImageStack(stack,100);

    for f = 1:size(stack,1)
        imwrite(stack{f},fullfile(path,strcat(name, '-', num2str(f), '.tif')),'tiff');
    end

end

