classdef AugmentedOptiSol < handle
% --------------------------------------------------------------------------
% AugmentedOptiSol
%   Wrapper class for casadi.OptiSol. 
%   Used as output of AugmentedOpti.solve_NLPSOL
% 
%
% METHODS:
% From casadi.OptiSol
%   - value -
%   - stats -
% 
%
% Original author: Lars D'Hondt
% Original date: 12/April/2023
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

    properties
        % Solution (casadi.OptiSol object)
        sol
        % Stats
        sol_stats
    end % end of public properties

    properties(Access = private)
        % Variables in cell array
        opti_vars
        % Variables in table
        opti_vars_table
    end % end of private properties

    methods
        % Constructor
        function obj = AugmentedOptiSol(sol,opti_vars,varargin)
            import casadi.*
            obj.sol = sol;
            obj.opti_vars = opti_vars;
            obj.opti_vars_table = struct2table(opti_vars);
            if ~isempty(varargin)
                obj.sol_stats = varargin{1};
            elseif isa(sol,'casadi.OptiSol')
                obj.sol_stats = sol.stats;
            end
        end

        % Read values
        function val = value(obj,opti_var)
            import casadi.*
            if isa(obj.sol,'casadi.OptiSol')
                val = obj.sol.value(opti_var);
            else
                val = value_double(obj,opti_var);
            end
        end
        
        % Get stats
        function sol_stats = stats(obj)
            sol_stats = obj.sol_stats;
        end
    end % end of public methods

    methods(Access = private)
        % Read values
        function val = value_double(obj,opti_var)
            import casadi.*

            opti_var_class = opti_var.class_name();

            if contains(opti_var_class,'Symbolic')
                % if defined as opti_var = opti.variable(...)
                opti_var_name = opti_var.name;
            else
                % if defined as opti_var = double*opti.variable(...)
                opti_var_name = opti_var.dep(1).name();
            end

            idx = find(strcmp(obj.opti_vars_table.name,opti_var_name));

            if isempty(idx)
                val = "not found";
            else
                val_vec = obj.sol(obj.opti_vars_table.index{idx});
                [val_size] = obj.opti_vars_table.size{idx};
                val = reshape(val_vec,val_size(1),val_size(2));
            end
        end
    end % end of private methods

end % end of class