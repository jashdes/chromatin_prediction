function [] = postdenoise(loc, template)
    template = imread(template);

    % Get list of files matching our naming scheme
    files = dir(fullfile(loc));
    % Remove non-directories from listing
    files = files([files.isdir]);
    files=files(~ismember({files.name},{'.','..', 'track'}));
    
    for f = 1 : size(files,1)
        % Get folder name
        foldername = files(f).name;
        
        care = fullfile(loc,foldername, 'dn');
        gtpath = fullfile(loc,foldername,strcat(foldername, '_gt.tif'));
        metric = isfile(gtpath);
        
        if metric
            gt = loadImageStack(gtpath);
            cropframe = imread(fullfile(loc,foldername,'gt',"cropframe.tif"));
        else
            cropframe = imread(fullfile(loc,foldername,"cropframe.tif"));
        end
        
        % Assemble stacks
        
        care = assembleStackDenoise(care);
        [~,~,~,care] = findGrid(cropframe,template,care);
        
        % Normalize
        for k = 1:size(care,1)
            care{k} = care{k} - min(min(care{k}));
            care{k} = care{k} ./ max(max(care{k}));
        end
        
        saveImageStack(care, fullfile(loc,foldername, strcat(foldername,'_dn.tif')));
        
        if metric
            [~,~,~,gt] = findGrid(cropframe,template,gt);
            saveImageStack(gt, fullfile(loc,foldername, strcat(foldername,'_gtc.tif')));
            % Copy all stacks to output folder
            copyfile(fullfile(loc,foldername,strcat(foldername,'_gtc.tif')), fullfile(loc,strcat(foldername,'_gt.tif')));
        end
        
        copyfile(fullfile(loc,foldername,strcat(foldername,'_dn.tif')), fullfile(loc,strcat(foldername,'_dn.tif')));
        
    end
end
