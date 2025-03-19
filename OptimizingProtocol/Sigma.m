classdef Sigma
    % Implement Pauli matrices i.e. SU2 algebra
    
    properties
        i
    end
    
    methods
        function obj = Sigma(ind)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.i = ind;
        end
        
        function str = r(obj)
            str = ['s',num2str(obj.i)];
        end
        
        
    end
end

