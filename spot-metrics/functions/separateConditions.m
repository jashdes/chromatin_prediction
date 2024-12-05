function [subsets] = separateConditions(d)
%SEPARATECONDITIONS Organizer to separate conditions in a directory listing

    %d:         directory listing of all files
    %subsets:   cell array with first column denoting output names
    %           and second column containing sub-listings
    
    % trim folder entries
    d = d([d.isdir] == 0);
    
    % initialize subsets
    map = containers.Map('KeyType','char','ValueType','any');
    unknown = struct();
    
    % condition notation convention
    expression = '_*cell\d+';
    
    % for each folder found
    for j = 1 : size(d,1)
        
        % Get name
        filename = d(j).name;
        
        % Extract condition string
        fparts = regexpi(filename,expression,'split');
        
        % If there is only one element, the naming convention has not been
        % followed
        if length(fparts) < 2
            warning("%s does not follow naming convention\n", filename);
            % append to unknown
            if isempty(fieldnames(unknown))
                unknown = d(j);
            else
                unknown = [unknown; d(j)];
            end
        else
            if isKey(map, fparts{1})
                list = map(fparts{1});
                list = [list; d(j)];
                map(fparts{1}) = list;
            else
                map(fparts{1}) = d(j);
            end
        end
    end
    
    % package subsets
    if ~isempty(fieldnames(unknown))
        subsets = {"unknown", unknown};
    else
        subsets = {};
    end
    
    % add all found conditions to subset list
    ks = keys(map);
    for k = 1 : length(ks)
        subsets(end+1,:) = {ks{k}, map(ks{k})};
    end

end

