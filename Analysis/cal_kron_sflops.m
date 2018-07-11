function [flops_total, flops_method1, flops_method2, flops_method3] = ...
      cal_kron_sflops( nrowA,nrowB, ncolA, ncolB, nnzA,nnzB )
% [flops_total, flops_method1, flops_method2, flops_method3]  = ...
%     cal_kron_sflops( nrowA,ncolA, nrowB, ncolB, nnzA,nnzB )
%
% calcuate the number of SPARSE flops needed to compute
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
%            (ii) Y += BX * transpose(A) or (A * transpose(BX))^t
% -------------------------------------- 
% flops_BX = 2*nrowB * ncolB * ncolX;
flops_BX = 2*nnzB * ncolX;

ncolBX = ncolX;
nrowBX = nrowB;
%flops_BX_At =  2*nrowY.*ncolY .* ncolBX;
flops_BX_At = 2*nnzA * nrowBX;

flops_method1 = flops_BX + flops_BX_At;

% ---------------------------------------------
% Method 2:  (i) compute XAt = X * transpose(A) or (A * transpose(X))^t
%            (ii) Y += B * XAt
% ---------------------------------------------
%flops_XAt = 2*nrowX .* ncolX .* nrowA;
flops_XAt =  2*nnzA * nrowX;

ncolXAt = nrowA;
% flops_B_XAt = 2 * nrowB .* ncolB .* ncolXAt;
flops_B_XAt = 2*nnzB * ncolXAt;

flops_method2 = flops_XAt + flops_B_XAt;

% ------------------------------
% Method 3: expand kron(A,B) * X
% ------------------------------
flops_method3 = 2 * nnzA * nnzB;

flops_total = nrowY.*ncolY + min( flops_method1, min(flops_method2, flops_method3) );

