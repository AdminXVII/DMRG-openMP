function [flops_total, flops_method1, flops_method2] = cal_kron_flops( nrowA,nrowB, ncolA, ncolB )
% [flops_total, flops_method1, flops_method2]  = cal_kron_flops( nrowA,ncolA, nrowB, ncolB )
% calcuate the number of flops needed to compute
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
flops_BX = 2*nrowB * ncolB * ncolX;

ncolBX = ncolX;
flops_BX_At =  2*nrowY.*ncolY .* ncolBX;
flops_method1 = flops_BX + flops_BX_At;

% ---------------------------------------------
% Method 2:  (i) compute XAt = X * transpose(A)
%            (ii) Y += B * XAt
% ---------------------------------------------
flops_XAt = 2*nrowX .* ncolX .* nrowA;

ncolXAt = nrowA;
flops_B_XAt = 2 * nrowB .* ncolB .* ncolXAt;
flops_method2 = flops_XAt + flops_B_XAt;

flops_total = nrowY.*ncolY + min( flops_method1, flops_method2 );

