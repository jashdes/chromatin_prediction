function [] = exportCellData(data, outpath)
%EXPORTDATA Exports batch data as specified in forEachFile for whole cells

%   data:       n x 2 cell where n is the number of processed files,
%               the first column is the file name, the second the results
%   outpath:    path to which results should be written

    % set default filename if none provided
    if(nargin < 2)
        outpath = 'results.csv';
    else
        % add .csv extension if not passed
        [~,~,ext] = fileparts(outpath);
        if ~strcmpi(ext,'.csv')
            outpath = strcat(outpath,'.csv');
        end
    end

    % open output file
    file = fopen(outpath,'w');
    
    k = 1;
    
    % get data struct
    filedata = data{:,2};

    % write file name to output file
    fprintf(file,'filename,');

    % get list of fields in the data struct
    fields = fieldnames(filedata);

    % write field names to first row
    for i=1:numel(fields)

      fieldname = fields{i};
      if i < numel(fields)
          fieldname = strcat(fieldname, ',');
      end
      fprintf(file,'%s', fieldname);

    end

    % for reach file entry
    for k = 1 : size(data,1)
        
        % get data struct
        filedata = data{k,2};
        
        % start next row
        fprintf(file,'\n');
        
        fprintf(file,'%s,', data{k,1});
        
        % for each row of the longest filedata entry
        for j = 1 : max(structfun(@numel,filedata))
            % for each filedata entry
            for i=1:numel(fields)
                % if no elements are left to be written
                if j > numel(filedata.(fields{i}))
                    % create empty dummy datum
                    datum = '';
                else
                    % read datum from field
                    datum = num2str(filedata.(fields{i})(j));
                end
                
                % append comma if not at last field
                if i < numel(fields)
                    datum = strcat(datum,',');
                end
                
                % write datum
                fprintf(file,"%s", datum);
            end
            
        end
        
    end
    
    % close file
    fclose(file);

end
