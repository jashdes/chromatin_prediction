function [] = saveImageStack(stack,path)
%SAVEIMAGESTACK Saves cell array representing image stack to .tiff

    %   stack:  cell array containing the image stack
    %   path:   path to which the .tiff will be written
    
    % for every frame
    for f = 1 : size(stack,1)
        
        % if this is the first frame
        if f == 1
            % create new file and write first frame
            imwrite(stack{f},path,'tiff');
        else
            % append current frame to existing stack
        imwrite(stack{f},path,'tiff','WriteMode','append');
        end
    end
            
    
end

