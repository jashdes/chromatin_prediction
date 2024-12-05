function [registered] = registerStackByOutline(unregistered)
%REGISTERSTACKBYOUTLINE Registers imagestack by applying StackReg to nucleus outline

    % segment image stack
    segmented = cellfun(@(img) segmentImage(img),unregistered,'un',0);
    
    % cast to int16 for StackReg
    segmented = cellfun(@im2uint16, segmented, 'un', 0);
    
    % register segmentations using StackReg
    registered = stackReg(segmented);
    
    % for every frame
    for frame = 1:length(registered)
        % extract transformation made by StackReg
        reg = registerImages(segmented{frame},registered{frame}, ...
            unregistered{frame});
        
        % and apply it to the unregistered image
        registered{frame} = reg.samewarp;
    end

end

