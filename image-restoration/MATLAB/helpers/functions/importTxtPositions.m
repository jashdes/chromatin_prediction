function [pos] = importTxtPositions(path)
%IMPORTTXTPOSITIONS Imports positions from Dr. Holzwarth's tracker output
% path:     path to .txt file containing positions as follows:
%           n lines per spot, where n is the number of dimensions
%           f+1 columns, where f is the number of frames being tracked
%           the first column contains the spot number, for the remainder,
%           column f+1 contains frame f's coordinates

file = readmatrix(path);

% drop first column containing spot numbers
file(:,1) = [];
nframes = size(file,2);
nspots = size(file,1) / 2;
assert (floor(nspots) == nspots);

% pre-allocate result array
pos = cell(nframes,1);
    
% for each frame
    
for f = 1:nframes

    % write results to frame struct
    frame.x = file(1:2:end,f); % odd rows
    frame.y = file(2:2:end,f); % even rows
    
    % write frame struct to result array

    pos{f} = frame;

end

end

