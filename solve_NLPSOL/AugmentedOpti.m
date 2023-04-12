classdef AugmentedOpti < handle
% --------------------------------------------------------------------------
% AugmentedOpti
%   Wrapper class for casadi.Opti. This class offers the possibility to use
%   optimisation variable bounds as bounds in the interior point optimisation. 
%   The regular optistack, with Opti.solve(), considers bounds as inequality 
%   constraints. AugmentedOpti uses the same syntax as Opti.
% 
%
% METHODS:
% Specific to AugmentedOpti
%   - AugmentedOpti() -
%   * Constructor. Returns an AugmentedOpi object
%
%   - solve_NLPSOL() -
%   * Solves the optimisation problem, but with bounds. Returns an
%   AugmentedOptiSol object.
%
% From casadi.Opti
%   - variable -
%   - parameter -
%   - set_initial -
%   - set_value -
%   - subject_to -
%   - minimize -
%   - solver -
%   - solve -
%   
%   To call other methods, use (obj).opti.(method)
% 
%
% Original author: Lars D'Hondt
% Original date: 12/April/2023
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

    properties
        % casadi.Opti object
        opti 
    end % end of public properties

    properties(Access = private)
        % Number of optimisation variables
        N_vars = 0;
        % Cell array with information about optimisation variables
        opti_vars = [];
        % Solver name
        solver_name
        % Solver options
        solver_options
    end % end of private properties

    methods
        % Constructor
        function obj = AugmentedOpti()
            import casadi.*
            obj.opti = casadi.Opti();
        end

        %% Problem formulation
        % Add a variable
        function opti_var = variable(obj,varargin)
            import casadi.*
            opti_var = obj.opti.variable(varargin{:});
            obj.opti_vars(end+1).name = opti_var.name();
            obj.opti_vars(end).var = opti_var;
            obj.opti_vars(end).size = {size(opti_var)};
            idx = [1:numel(opti_var)] + obj.N_vars;
            obj.N_vars = idx(end);
            obj.opti_vars(end).index = idx;
        end

        % Add a parameter
        function [] = parameter(obj,varargin)
            import casadi.*
            obj.opti.parameter(varargin{:});
        end

        % Set initial guess
        function [] = set_initial(obj,varargin)
            import casadi.*
            obj.opti.set_initial(varargin{:});
        end

        % Set parameter value
        function [] = set_value(obj,varargin)
            import casadi.*
            obj.opti.set_value(varargin{:});
        end

        % Add constraint
        function [] = subject_to(obj,varargin)
            import casadi.*
            obj.opti.subject_to(varargin{:});
        end

        % Set objective
        function [] = minimize(obj,varargin)
            import casadi.*
            obj.opti.minimize(varargin{:});
        end

        %% Solve
        % Set solver
        function [] = solver(obj,name,options)
            import casadi.*
            obj.opti.solver(name,options);
            obj.solver_options = options;
            obj.solver_name = name;
        end

        % Call solver (opti.solve) that considers bounds as inequality
        % contraints
        function [sol,varargout] = solve(obj)
            import casadi.*
            sol = obj.opti.solve();
            if nargout == 2
                varargout{1} = AugmentedOptiSol(sol,obj.opti_vars);
            end
        end

        % Call solver that considers bounds as bounds
        function [aug_sol] = solve_NLPSOL(obj)
            import casadi.*

            [sol,stats] = solve_NLPSOL_for_AugmentedOpti(obj.opti,...
                obj.solver_options, obj.solver_name);

            aug_sol = AugmentedOptiSol(sol,obj.opti_vars,stats);

        end

        %% Get private properties
        % Number of optimisation variables
        function N_vars = get_total_number_of_variables(obj)
            N_vars = obj.N_vars;
        end

        % Cell array with information about optimisation variables
        function opti_vars = get_variables_bookkeeping(obj)
            opti_vars = obj.opti_vars;
        end

        % Solver name
        function solver_name = get_solver_plugin_name(obj)
            solver_name = obj.solver_name;
        end

        % Solver options
        function solver_options = get_solver_options(obj)
            solver_options = obj.solver_options;
        end
    end % end of public methods
    
end % end of class