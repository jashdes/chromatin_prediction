function [] = exportDenoise(data, outpath, fixed)
%EXPORTDATA Exports batch data as specified in forEachFile for whole cells

%   data:       path to experiment
%   outpath:    path to which results should be written

    % set default filename if none provided
    if(nargin < 2)
        outpath = data;
    end
    
    if(nargin < 3)
        fixed = false;
    end
    
    nspots = 49;
    cells = getDirs(fullfile(data,'data'));
    methods = getDirs(fullfile(data,'track'));
    celltable = [];
    spottable = [];
    
    for m = 1:length(methods)
        mname = methods(m).name;
        [~, res] = stackRestorationMetrics(fullfile(data,'track',mname,'txt'), fixed);
    
        for c = 1:length(cells)
            % per-cell
            cname = cells(c).name;
            cres = res(strcmp({res.name}, cname));
            cSSdnpath=fullfile(data,'track',mname,'mat',strcat(cname,"_dn_SS.mat"));
            cSSgtpath=fullfile(data,'track',mname,'mat',strcat(cname,"_gt_SS.mat"));
            wasTracked = isfile(cSSdnpath) && (isfile(cSSgtpath) || fixed);
            if wasTracked
                cSS = load(cSSdnpath);
                cSS = cSS.SS;
                dmean = cSS.field_Dmean;
                spotstracked = cres.nspots;
                posdeltas = cres.posdeltas;
                
            else
                dmean = NaN;
                spotstracked = NaN;
                posdeltas = NaN;
                
            end
            ctable = table(string(cname),string(mname), dmean, spotstracked, posdeltas, ...
                    'VariableNames',["cell","method","Dmean","nspots", "trackerror"]);
            if c == 1 && m == 1
                celltable = ctable;
            else
                celltable = [celltable; ctable];
            end
                

            % per-spot
            spotn = transpose(1:nspots);
            if wasTracked
                D = cSS.field_D;
                dD = cSS.field_dD;
            else
                D = NaN(nspots,1);
                dD = NaN(nspots,1);
            end
            
            stable = table(repmat(string(cname),nspots,1),repmat(string(mname),nspots,1), spotn, D, dD, ...
                    'VariableNames',["cell","method","spotN","D", "dD"]);
            
            if c == 1 && m == 1
                spottable = stable;
            else
                spottable = [spottable; stable];
            end

        end
    
    end
    
    [~,expname,~] = fileparts(data);
    writetable(celltable,fullfile(outpath, strcat(expname,"_cells.csv")));
    writetable(spottable,fullfile(outpath, strcat(expname,"_spots.csv")));

end
