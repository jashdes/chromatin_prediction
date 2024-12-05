function [found,missing] = separateFilesByName(d,s)
%SEPARATEFILESBYNAME Organizes directory listing into files that include a
% substring and those that do not

    %   d:  directory listing
    %   s:  substring to determine association
    
    indices = contains({d.name},s,'IgnoreCase',true);
    found = d(indices);
    missing = d(~indices);

end

