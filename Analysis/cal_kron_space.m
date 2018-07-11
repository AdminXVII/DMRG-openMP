function [space_method1, space_method2] = cal_kron_space( nrowA,nrowB, ncolA, ncolB )
% [space_total, space_method1, space_method2]  = cal_kron_space( nrowA,ncolA, nrowB, ncolB )
% calcuate the amount of temporary space needed to compute
%
% Y += kron(A,B) * X 
% or
% Y += B * X * transpose(A)
%
% where matrix A is nrowA by ncolA
% and matrix B is nrowB by ncolB
%

nrowX = ncolB;
ncolX = ncolA;

nrowY = nrowB;
ncolY = nrowA;
% -------------------------------------- 
% Method 1:  (i)  compute BX = B*X, 
%            (ii) Y += BX * transpose(A)
% -------------------------------------- 
space_BX = nrowB * ncolX;
space_method1 = space_BX;


% ---------------------------------------------
% Method 2:  (i) compute XAt = X * transpose(A)
%            (ii) Y += B * XAt
% ---------------------------------------------
space_XAt = nrowX * nrowA;
space_method2 = space_XAt;


