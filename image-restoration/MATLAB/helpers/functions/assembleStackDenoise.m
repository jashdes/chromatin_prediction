function [stack] = assembleStackDenoise(path)
%ASSEMBLESTACKDENOISE Re-assembles denoised exploded image stack

%   path:   folder containing individual denoised images

    stack = cell(1,1);

    % Get list of files matching our naming scheme
    files = dir(fullfile(path, '*-*.tif'));
    % Remove directories from listing
    files = files(~[files.isdir]);
    
    for f = 1 : size(files,1)
        % Get file name
        filename = files(f).name;
            
        % Extract frame name
        frame = regexpi(filename,'-(\d*).tif$', 'tokens');
        if all(size(frame) == [1,1])
            frame = str2num(frame{1}{1});
        else
            error("Frame number not found");
        end
        
        stack{frame,1} = im2double(imread(fullfile(path, filename)));
        
    end
    
end

