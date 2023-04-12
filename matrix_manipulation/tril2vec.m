function [vec] = tril2vec(tril)
% --------------------------------------------------------------------------
% tril2vec
%   Write the non-zero elements of a square lower triangular matrix into a
%   column vector.
%   To evaluate this function for multiple inputs, concatenate them
%   horizontally in a matrix.
%
%   Example:
%       tril = [a  0  0
%               d  b  0
%               e  f  c];
%       tril = vec2tril(vec)
%       vec = [a b c d e f]'
%       
% Reverso operation: vec2tril
%
% INPUT:
%   - tril -
%   * square, lower triangular matrix.
%
% OUTPUT:
%   - vec -
%   * a column vector
% 
%
% Original author: Lars D'Hondt
% Original date: 07/November/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

% Test input size
n = size(tril,1);
m = size(tril,2);
if rem(m,n) ~= 0
    error('Expected a square matrix as input.')
end

% Define vector
N = n*(n+1)/2;

% Call implementation
if n==m
    [vec] = tril2vec_impl(tril,n,N);
else
    [vec] = tril2vec_horzcat(tril,n,N,m);
end

% If result is a casadi double, return as a matlab double
if isa(vec,'casadi.DM')
    vec = full(vec);
end

end % end of function tril2vec

%% Implementation of tril2vec for column vector
function [vec] = tril2vec_impl(tril,n,N)
    import casadi.*
%     % decide variable type based on input
%     if isa(tril,'casadi.SX')
%         vec = SX(N,1);
%     elseif isa(tril,'casadi.MX')
%         vec = MX(N,1);
%     elseif isa(tril,'casadi.DM')
%         vec = DM(N,1); 
%     else
%         vec = zeros(N,1);
%     end

    vec = repmat(tril(1,end),N,1);

    % Get diagonal elements
    vec(1:n,1) = diag(tril);
    i_vec = n;
    
    % Get elements below diagonal
    for i=1:n-1
        idx = i_vec + [1:i];
        vec(idx,1) = tril(i+1,1:i);
        i_vec = idx(end);
    end

end % end of function tril2vec_impl

%% Repeat tril2vec for each tril of vec, and horizontally concatenate columns of vec.
function [vec] = tril2vec_horzcat(tril,n,N,m)
    import casadi.*
    % create casadi function of tril2vec
    tril_SX = SX.sym('tril_SX',Sparsity.lower(n));
    vec_SX = tril2vec_impl(tril_SX,n,N);
    f_tril2vec = Function('tril2vec',{tril_SX},{vec_SX},{'tril_SX'},{'vec_SX'});
    % map function to number of columns in vec (=m)
    f_tril2vec_concat = f_tril2vec.map(m/n);
    % evaluate mapped function to get all columns of vec
    vec = f_tril2vec_concat(tril);

end % end of function tril2vec_horzcat