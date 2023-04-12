function [tril] = vec2tril(vec)
% --------------------------------------------------------------------------
% vec2tril
%   Reshape a vector to a square lower triangular matrix. 
%   To evaluate this function for multiple inputs, concatenate them
%   horizontally in a matrix.  
%
%   Example:
%       vec = [a b c d e f]';
%       tril = vec2tril(vec)
%       tril = [a  0  0
%               d  b  0
%               e  f  c]
% 
% Reverse operation: tril2vec
%
% INPUT:
%   - vec -
%   * a column vector
%
% OUTPUT:
%   - tril -
%   * square, lower triangular matrix.
% 
%
% Original author: Lars D'Hondt
% Original date: 07/November/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

% Find size for nxn matrix
N = size(vec,1);
% Solve N = n*(n+1)/2 for n
n = max(roots([1/2,1/2,-N]));

% Test size
if rem(n,1) ~= 0
    n_next = ceil(n);
    N_needed = n_next*(n_next+1)/2;
    error(['A lower triangular ' num2str(n_next) 'x' num2str(n_next) ' matrix has ',...
        num2str(N_needed) ' non-zero elements. The provided vector has ' num2str(N) ' elements.'])
end

% Call implementation
if size(vec,2) == 1
    [tril] = vec2tril_impl(vec,n);
else
    [tril] = vec2tril_horzcat(vec,n,N);
end

% If result is a casadi double, return as a matlab double
if isa(tril,'casadi.DM')
    tril = full(tril);
end

end % end of function vec2tril

%% Implementation of vec2tril for column vector
function [tril] = vec2tril_impl(vec,n)
    % Define square matrix and set diagonal elements
    tril = diag(vec(1:n,1));
    i_vec = n;
    
    % Set elements below diagonal;
    for i=1:n-1
        idx = i_vec + [1:i];
        tril(i+1,1:i) = vec(idx,1);
        i_vec = idx(end);
    end
end % end of function vec2tril_impl

%% Repeat vec2tril for each column of vec, and horizontally concatenate the tril.
function [tril] = vec2tril_horzcat(vec,n,N)
    import casadi.*
    % create casadi function of vec2tril
    vec_SX = SX.sym('vec_SX',N,1);
    tril_SX = vec2tril_impl(vec_SX,n);
    f_vec2tril = Function('vec2tril',{vec_SX},{tril_SX},{'vec_SX'},{'tril_SX'});
    % map function to number of columns in vec
    f_vec2tril_concat = f_vec2tril.map(size(vec,2));
    % evaluate mapped function to get all columns of vec
    tril = f_vec2tril_concat(vec);

end % end of function vec2tril_horzcat