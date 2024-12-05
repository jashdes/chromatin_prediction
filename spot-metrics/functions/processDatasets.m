function [err] = processDatasets(path,analyses, outpath)
%PROCESSDATASETS Performs analyses on datasets in specified path

    %   path:           path pointing to directory that contains datasets
    %   analyses:       array of analysis structs
    %   analysis
    %   .name           descriptive name of analysis
    %   .foldername:    name of folder containing relevant data
    %   .organizer:     function to organize data before passing to handler
    %   .handler:       function to be called to do this analysis with a
    %                   path to each file 
    %   .args:          arguments passed to handler
    %   outpath:        path for the results to be written to
    
    % get a directory listing
    
    d = dir(path);
    
    % trim file rows
    d = d([d.isdir] == 1);
    
    % for each folder found
    for j = 1 : size(d,1)
        
        % skip UNIX directory symbols
        if strcmp(d(j).name, '.') || strcmp(d(j).name, '..')
            continue;
        end
        
        %   for each analysis
        for k = 1 : size(analyses,1)

            analysis = analyses(k);

            % determine full paths
            
            % for required folder
            p = fullfile(path,d(j).name,analysis.foldername);
            
            % for results
            resp = fullfile(outpath,d(j).name,analysis.name);

            % check if the required folder exists
            if exist(p,'dir')
                
                % get a directory listing to pass
                adir = dir(p);
                
                % trim folder rows
                adir = adir([adir.isdir] == 0);
                
                % create output directory
                mkdir(resp);
                
                % apply organizer if passed
                if isfield(analysis,'organizer')
                    subsets = analysis.organizer(adir);
                    
                    % perform analysis for each subset
                    for s = 1 : size(subsets,1)
                        % generate output path
                        op = fullfile(resp,subsets{s,1});
                        
                        forEachFile(subsets{s,2}, analysis.handler, ...
                            analysis.args, op);
                    end
                else
                    % perform analysis on complete directory listing
                    forEachFile(adir,analysis.handler,analysis.args, resp);
                end
                
            else
                warning('No %s directory; skipping %s...', ...
                    analysis.foldername, ...
                    analysis.name);
            end
        end
    end
    
    
    
end
