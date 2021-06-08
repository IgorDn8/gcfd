function A = vec_mismatch (V1, V2)
% Vector Mismatch Measure 
%
% Takes values from 0 to 1. 
% 
% INPUT: 
%     V1, V2 - two matrices of the same size. 
%
% OUTPUT:
%     A -     the mismatch measure. For collinear vectors, A == 0. For
%     orthogonal vectors, A == 1. The operation is executed column-wise

A = 1-abs(dot(V1,V2,1))./(sqrt(sum(V1.^2,1)).*sqrt(sum(V2.^2,1)));