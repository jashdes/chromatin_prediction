function [subsets] = separateControls(d)
%SEPARATECONTROLS Organizer to separate controls from other files

    %d:         directory listing
    %subsets:   cell array with first column denoting output names
    %           and second column containing sub-listings
    
    [control,experimental] = separateFilesByName(d,'control');
    
    subsets = {'control',control;'experimental',experimental};

end

