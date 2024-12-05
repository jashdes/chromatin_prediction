function [registered] = stackReg(unregistered)
%STACKREG Applies ImageJ StackReg on image set

    %   unregsitered:   unregistered images tack
    %   registered:     result of applying StackReg to unregistered
    
    % Requires http://bigwww.epfl.ch/sage/soft/mij/mij.jar in javapath
    % Modified from https://www.mathworks.com/matlabcentral/fileexchange/45633-stackregwrapper
    
    % Prepare a headless ImageJ
    Miji(false);
    
    % remember data type
    dtype = class(unregistered{1});
    
    % Convert cell array image stack to 3D matrix for transfer
    unregistered = cat(3,unregistered{:});
    
    % Transfer image stack to ImageJ
    MIJ.createImage(unregistered);
   
    % Register
    MIJ.run("StackReg ", "transformation=[Rigid Body]");
    
    % transfer image from FIJI to Matlab
    registered = MIJ.getCurrentImage;
    
    % close image window in FIJI
    MIJ.run('Close');

    % exit FIJI
    MIJ.exit
    
    % convert image stack back to cell array
    registered = num2cell(registered, [1,2]);
    
    % remove zero dimensions
    registered = squeeze(registered);
    
    % cast image stack to original data type
    registered = cellfun(str2func(dtype),registered, 'un', 0);
    
end

