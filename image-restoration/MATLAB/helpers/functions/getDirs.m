function [folders] = getDirs(loc)

    % Get list of files matching our naming scheme
    folders = dir(loc);

    % Remove non-directories from listing
    folders = folders([folders.isdir]);
    folders = folders(~ismember({folders.name},{'.','..'}));
    
end

