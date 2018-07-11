#include "analysis.h"
void cal_kron_flops( 
        int nrowA, int nrowB,
        int ncolA, int ncolB,
        double *pflops_total,
        double *pflops_method1,
        double *pflops_method2 )
{

 int nrowX = ncolB;
 int ncolX = ncolA;

 int nrowY = nrowB;
 int ncolY = nrowA;


    /*
    % -------------------------------------- 
    % Method 1:  (i)  compute BX = B*X, 
    %            (ii) Y += BX * transpose(A)
    % -------------------------------------- 
    */
    double flops_BX = (2.0*nrowB) * ncolB * ncolX;
    
    int ncolBX = ncolX;
    double flops_BX_At =  (2.0*nrowY) * ncolY * ncolBX;
    double flops_method1 = flops_BX + flops_BX_At;
    
    /*
    % ---------------------------------------------
    % Method 2:  (i) compute XAt = X * transpose(A)
    %            (ii) Y += B * XAt
    % ---------------------------------------------
    */
    double flops_XAt = (2.0*nrowX)  * ncolX  * nrowA;
    
    int ncolXAt = nrowA;
    double flops_B_XAt = (2.0 * nrowB)  * ncolB  * ncolXAt;
    double flops_method2 = flops_XAt + flops_B_XAt;
    
    double flops_total = nrowY *ncolY + MIN( flops_method1, flops_method2 );

    *pflops_total = flops_total;
    *pflops_method1 = flops_method1;
    *pflops_method2 = flops_method2;
}

