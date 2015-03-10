classdef murates
    % Motor unit firing rate curves
    
    properties
        RET
        Rmin
        Rmax
        m
    end
    
    methods
        function obj = murates(RET,Rmin,Rmax,m)
            obj.RET = RET;
            obj.Rmin = Rmin;
            obj.Rmax = Rmax;
            obj.m = m;
        end
        
        function rate = getRate(obj,exc,N)
            deltaexc = exc - obj.RET(N);
            if deltaexc < 0
                rate = 0;
            else
                rate = min(obj.Rmin + obj.m*deltaexc,obj.Rmax(N));
            end
        end
    end
end