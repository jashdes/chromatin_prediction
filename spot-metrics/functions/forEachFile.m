function [results, errors] = forEachFile(d, f, args, outpath, wholecell, shouldpause)
%FOREACHFILE Runs function f for each file in folder path

    % d:            directory listing as obtained from calling dir(path)
    % f:            function to process the files; will receive absolute path of
    %               each file as argument
    % args:         argument to be passed after the path when f is called
    % args.wholecell:set if the output file summarizes cell properties
    % args.writer:  set if f contains a function to write results
    %               internally
    % outpath:      path to which results will be written if passed
    % shouldpause:  will pause after each execution if set to 1
    % results:      cell array of file names and structures as returned by f
    % errors:       cell array of files names and MExceptions as thrown by calls to f
    
    % trim folder rows
    d = d([d.isdir] == 0);
    
    if(nargin < 6)
        shouldpause = 0;
    end
    
    % pre-allocate result cell
    results = cell(size(d,1),2);
    
    % for each file
    for k = 1:size(d,1)
        % report progress
        fprintf("%4.2f%%: processing %s\n", ((k-1)/size(d,1))*100,d(k).name);

        % get the absolute path to the file
        filepath = fullfile(d(k).folder,d(k).name);
        
        try
            % write file name to results
            results{k,1} = d(k).name;
            
            % if no arguments have been given to parse
            if(nargin < 3)
                % call f with its absolute path as the sole argument
                results{k,2} = f(filepath);
            else
                if isfield(args,'writer')
                    % add output path to arguments to be used by writer
                    args.outpath = outpath;
                end
                % else call f with the absolute path and arguments
                results{k,2} = f(filepath, args);
            end
            
        catch err
            warning('An error occured while processing %s\n%s\n', ...
                d(k).name, err.message);
            results{k,2} = err;
            
        end
        if(shouldpause)
            pause on;
            pause;
        end
    end
    
    % separate errors and results
    err = cell(size(results,1),1);
    err(:) = {'MException'};
    iserror = cellfun(@strcmp, ...
        cellfun(@class,results(:,2), 'UniformOutput', false),err);
    errors = results(iserror,:);
    results = results(~iserror,:);
    
    % rely on writer if provided, otherwise write to .csv
    if(nargin >= 4 && ~isfield(args,'writer'))
        if(isfield(args,'wholecell') && args.wholecell)
            exportCellData(results,outpath);
        else
            exportData(results, outpath);
        end
    end 
    
end
