function [total, res] = stackRestorationMetrics(d, fixed)
%STACKRESTORATIONMETRICS Measures discrepancy between gt and denoised stack

    if nargin < 1
        d = uigetdir('.','Pick folder containing .txt tracking data');
    end
    
    if nargin < 2
        fixed = false
    end

    ext = '_SS.mat.txt';
    res = {};
    total = 0;
    
    % Get list of files matching our naming scheme
    files = dir(fullfile(d, '*.txt'));
    % Remove directories from listing
    files = files(~[files.isdir]);
    cnt = 1;

    for f = 1 : size(files,1)
        % Get file name
        filename = files(f).name;

        % Check get base filename
        basename = regexpi(filename,strcat('(.*)_gt',ext), 'tokens');
        
        % Skip non ground truth tracks
        if size(basename,1) == 0
            continue
        end
        
        % Strip cell array
        basename = basename{1}{1};
        
        stack = struct();
    
        % load images from respective folders
        groundtruth = fullfile(d,strcat(basename, '_gt', ext));
        care = fullfile(d,strcat(basename, '_dn', ext));
        stack.imageset = {groundtruth; care};

        % test for number of trackable spots
        stack.trackableSpots = [0, 0];
        stack.positions = {NaN, NaN};
        stack.posdist = [NaN];

        for j = 1:2
            if ~isfile(stack.imageset{j})
                continue
            end
            stack.positions{j} = importTxtPositions( stack.imageset{j});
        end
        
        nframes = size(stack.positions{1},1);
        
        for frame = 1:nframes
            for j = 1:2
                
                if all(size(stack.positions{j}) == [1 1]) && isnan(stack.positions{j})
                    trackableSpots = 0;
                    continue;
                else
                    trackableSpots = size(stack.positions{j}{f}.x,1) - sum(isnan(stack.positions{j}{f}.x));
                end
                
                if frame == 1
                    stack.trackableSpots(j) = trackableSpots;
                else
                    stack.trackableSpots(j) = min(trackableSpots, stack.trackableSpots(j));
                end

                % compute tracking errors
                if j > 1
                    c = j - 1;
                    if frame == 1
                        stack.posdist(c) = 0;
                    else
                    xs = (stack.positions{1}{frame}.x-stack.positions{1}{frame-1}.x);
                    ys = (stack.positions{1}{frame}.y-stack.positions{1}{frame-1}.y);
                    
                    xd = (stack.positions{2}{frame}.x-stack.positions{2}{frame-1}.x);
                    yd = (stack.positions{2}{frame}.y-stack.positions{2}{frame-1}.y);
                    
                    xdiff = xs-xd;
                    ydiff = ys-yd;
                    
                        if ~fixed
                            stack.posdist(c) = stack.posdist(c) + sum(sqrt(xdiff.^2+ydiff.^2)/min(stack.trackableSpots(1),stack.trackableSpots(2)),'all','omitnan');
                        else
                            stack.posdist(c) = stack.posdist(c) + sum(sqrt(xd.^2+yd.^2)/min(stack.trackableSpots(1),stack.trackableSpots(2)),'all','omitnan');
                        end
                    end
                end

            end
        end

        res(cnt).name = basename;
        res(cnt).nspots = stack.trackableSpots(2)/stack.trackableSpots(1);
        res(cnt).posdeltas = stack.posdist/nframes;
        total = total + stack.posdist/nframes;
        cnt = cnt + 1;
        
    end
    nfiles = size(files,1);
    total = total/nfiles;
    
end