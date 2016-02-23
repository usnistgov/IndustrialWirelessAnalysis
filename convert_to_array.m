function convert_to_array(pattern)
% % Convert files to struct of meta + array of CIR's
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

more off; 
grid minor;
set(0,'DefaultFigureWindowStyle','docked')

%
cir_file = []; %#ok<NASGU>


files = dir(pattern);
for fk = 1:length(files)
    
    mat_fname = files(fk).name;
    disp(['opening ' mat_fname]);
    try 
        cir_file = load(mat_fname);
    catch me
        warning('Problem reading mat file, trying again then skipping.');
        disp(me.message)
        try
        cir_file = load(mat_fname);
        catch
            disp(me.message)            
            warning('Skipping file...');
            continue;
        end
    end
    
    the_run = struct('meta',[],'cir',[]);
    try
        the_run.meta = cir_file.AllthePdps(1:19,7:8);
    catch me
        warning('file appears empty. skipping...')
        disp(me.message);
        continue;
    end
    apf = the_run.meta{11,2};
    rpa = cell2mat(the_run.meta(9,2));
    ns =  cell2mat(the_run.meta(8,2));
    NN = apf*rpa;

    % setup output directories
    arr_dir = 'arr';
    mkdir(arr_dir);

    % form the CIR matrix
    the_run.cir = nan(ns,NN);
    for kk = 1:NN
        % extract the cir from cell array
        cir = cell2mat(cir_file.AllthePdps(kk,6));
        % save as column vectors
        if ~isempty(cir)
            the_run.cir(:,kk) = cir(:);
        end
    end
    save([arr_dir '\' mat_fname(1:end-4) '__array.mat'], 'the_run')

    % explicit clear of large memory
    cir_file = [];
    cir = [];
    
    % gather memory stats
    disp('memory usage after pack')
    memory

end

end % function
