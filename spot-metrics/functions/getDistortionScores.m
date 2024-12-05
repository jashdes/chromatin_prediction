function [res] = getDistortionScores(pos)
%GETDISTORTIONSCORES Get distortion scores for groups of tracking outputs
%   pos:    cell array of tracking outputs, shaped {filename, tracking}

    for k = 1:size(pos,1)
        [fLen, fAng] = cumulativeDistortion(0,0,pos{k,2}([1,end]));
        res{k,1} = pos{k,1};
        res{k,2}.ang = fLen(end);
        res{k,2}.len = fAng(end);
    end

end