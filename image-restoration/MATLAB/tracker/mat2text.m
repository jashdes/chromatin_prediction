function [] = mat2text(directory, out)
    if nargin == 0
        directory = uigetdir('.','Select folder containing .mat files');
        out = uigetdir('.','Select output folder');
    end

    % Get list of files matching our naming scheme
    files = dir(fullfile(directory, '*_SS.mat'));
    % Remove directories from listing
    files = files(~[files.isdir]);

    for k = 1 : size(files,1)
        % Get file name
        filename = files(k).name;
        pixelsPerMicron = 9.26;
        file = fullfile(directory,filename);

        load(file)
        A = SS.field_track_out;
        xROI_pixel = A(:,3);
        yROI_pixel = A(:,4);
        xROI_micron = A(:,3)/pixelsPerMicron;
        yROI_micron = A(:,4)/pixelsPerMicron;
        xAbsolutePixel = A(:,5);
        yAbsolutePixel = A(:,6);
        x = xAbsolutePixel;
        y = yAbsolutePixel;

        % Write data to file
        spotCount = 49;
        numFramesN = numel(x)/spotCount;
        M = zeros(numFramesN+1,2*spotCount);
        forCounter = 1;
        for i = 1:2:2*spotCount
            M(1,i:i+1) = (i+1)/2;
            M(2:numFramesN+1,i) = x(forCounter:forCounter+numFramesN-1);
            TF = isnan(M(2:numFramesN+1,i));
            if sum(TF) > 0
                M(2:numFramesN+1,i)=0; % set all values = 0 if even one original value is NaN
            end
            M(2:numFramesN+1,i+1) = y(forCounter:forCounter+numFramesN-1);
            TF = isnan(M(2:numFramesN+1,i+1));
            if sum(TF) > 0
                M(2:numFramesN+1,i+1)=0; % set all values = 0 if even one original value is NaN
            end
            forCounter = forCounter + numFramesN;
        end
        M(M==0) = NaN;
        file = fullfile(out, strcat(filename, ".txt"));
        dlmwrite(file,M','delimiter',' ');
    end
end
