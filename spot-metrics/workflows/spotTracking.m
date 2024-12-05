function [pos] = spotTracking(imstack, template, register)
%SPOTTRACKING Track spots within image sequences

    % imstack:  image stack representing sequence
    % pos:      struct cell array with fields x, y , sigmaX, sigmaY
    % register: whether the stack should be registered
    
    % pre-allocate result array
    pos = cell(size(imstack,1),1);
    
    % crop stack based on first image
    [~,imstack] = cropFullFrame(imstack{1},imstack);
    
    % register image stack if flag has not been passed or is set
    if register || nargin < 3
        imstack = stackReg(imstack);
    end
    
    % for each frame
    
    for f = 1:size(imstack,1)
        
        fprintf("%d/%d\n",f,size(imstack,1));
    
        % find the template grid within the image
        [grid,xoffset,yoffset] = findGrid(imstack{f},template);

        % find spot positions within the grid
        [gridpos, sigma] = findSpots(imadjust(grid),template);

        % calculate image positions from grid-relative positions
        framepos = gridpos + [yoffset,xoffset];
        
        % write results to frame struct
        frame.x = framepos(:,1);
        frame.y = framepos(:,2);
        frame.sigmaX = sigma(:,1);
        frame.sigmaY = sigma(:,2);
        
        % write frame struct to result array
        
        pos{f} = frame;
        
    end

end

