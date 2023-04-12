function [PD] = tril2PosDef(L)
% --------------------------------------------------------------------------
% tril2PosDef
%   Calculate a positive definite matrix PD based on a lower triangular
%   matrix L: PD = L*L'.
%   This functions uses the information that PD is symmetric to improve
%   efficiency.
%   To evaluate this function for multiple inputs, concatenate them
%   horizontally in a matrix.
% 
% Reverse operation: chol(PD,'lower')
%
% INPUT:
%   - tril -
%   * square, lower triangular matrix.
%
% OUTPUT:
%   - PD -
%   * positive definit matrix
% 
%
% Original author: Lars D'Hondt
% Original date: 07/November/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

import casadi.*

% Test input is square matrix
n = size(L,1);
m = size(L,2);
if rem(m,n) ~= 0
    error('Expected a square matrix as input.')
end

% Call implementation
if n==m
    PD = tril2PosDef_impl(L,n);
else
    PD = tril2PosDef_horzcat(L,n,m);
end

% If result is a casadi double, return as a matlab double
if isa(PD,'casadi.DM')
    PD = full(PD);
end

end % end of function tril2PosDef

%% Implementation
function [PD] = tril2PosDef_impl(L,n)
    import casadi.*
    % Cast input to sparse lower triangular matrix
    tril = project(L,Sparsity.lower(n));
    % Calculate PD as L*L'
    PD_tmp_1 = tril*tril';
    % Remove entries above diagonal
    PD_tmp_2 = project(PD_tmp_1,Sparsity.lower(n));
    % Add symmetry
    PD = tril2symm(PD_tmp_2);

end % end of function tril2PosDef_impl

%% Repeat tril2PosDef 
function [PD] = tril2PosDef_horzcat(L,n,m)
    import casadi.*
    % create casadi function of vec2tril
    L_SX = SX.sym('L_SX',Sparsity.lower(n));
    PD_SX = tril2PosDef_impl(L_SX,n);
    f_tril2PosDef = Function('vec2tril',{L_SX},{PD_SX},{'L_SX'},{'PD_SX'});
    % map function to number of columns in vec
    f_tril2PosDef_concat = f_tril2PosDef.map(m/n);
    % evaluate mapped function to get all columns of vec
    PD = f_tril2PosDef_concat(L);

end % end of function vec2tril_horzcat