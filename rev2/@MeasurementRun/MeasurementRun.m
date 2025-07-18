classdef MeasurementRun < handle
    
    properties(Constant)
        OPT_PATH_GAIN = 1;
        OPT_KFACTOR = 2;
        OPT_DELAY_SPREAD = 3;
        OPT_AVGCIR = 4;
        OPT_NTAP_APPROX = 5;
        OPT_WRITE_STATS = 6;
        OPT_DO_PLOTS = 7;     

        NtapApprox_N = 48;

        metrics_tbl_colnames = { ...
            'CoordX', 'CoordY', 'LOS', 'RicianK', ...
            'RMSDelaySpread', 'MeanDelay', 'MaxDelay', 'PathGain'};   
        metrics_tbl_vartypes = {'double', 'double', 'double', 'double', ...
            'double', 'double', 'double', 'double' };        

        % order of metrics
        CoordX = 1;
        CoordY =2;
        LOS =3;
        RicianK = 4;
        RMSDelaySpread = 5;
        MeanDelay = 6;
        MaxDelay = 7;
        PathGain = 8;
        % AverageCIR = 9:9+48-1; % needs to use NtapApprox_N
    end
    
    properties

        OPTS = zeros(7,1);

    end
    
    methods
        function obj = MeasurementRun(OPTS)
            if nargin < 1
                obj.OPTS = [ 1; 1; 1; 1; 1; 0; 0];
            else
                obj.OPTS = OPTS;
            end
        end
    end
end

