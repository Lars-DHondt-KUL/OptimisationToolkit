clear
close all
clc

import casadi.*

%% Test vec2tril and tril2vec

n = 3;
N = n*(n+1)/2;

% vec = SX.sym('vec',N,2);
% vec = MX.sym('vec',N,2);
syms a b c d e f g h real
vec = [a,b,c,d,e,f]';
vec_w = [g,h]';

tril = vec2tril(vec)

vec2 = tril2vec(tril);

diff = vec2 - vec

% P = tril2PosDef(tril)

P = tril*tril'
P1 = [P,zeros(n,2); zeros(2,n),diag(vec_w.^2)];

tril2 = [tril,zeros(n,2); zeros(2,n),diag(vec_w)];
P2 = tril2*tril2';

