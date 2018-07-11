#include <math.h>
#include <assert.h>
#include <stdio.h>

#include "analysis.h"

const int idebug = 0;

static double nchoosek(int n, int k) {
   k = MIN(k, n-k);

   if (k == 0) {
     return ( (double) 1.0);
     };
   
   double dval = 1.0;
   int i = 0;

   for(i=1; i <= k; i++) {
     double num = (double) (n-i+1);
     double den = (double) i;
     dval *= (num / den );
     };

   return (dval);
}

static double bincoeff(int n, int k) {
   return( nchoosek( n, k ) );
}



int gen_patches_comb( 
          int left_size,
          int right_size,
          int target_up, 
          int target_down, 
          int keep_left_states,
          int keep_right_states,

          int left_patch_size_[],
          int right_patch_size_[],

          int left_patch_up_[],
          int left_patch_down_[],
          int right_patch_up_[],
          int right_patch_down_[],

          int interaction_matrix_[]
           )
{
#define left_patch_size(i) left_patch_size_[(i)-1]
#define right_patch_size(i) right_patch_size_[(i)-1]

#define left_patch_up(i) left_patch_up_[(i)-1]
#define left_patch_down(i) left_patch_down_[(i)-1]

#define right_patch_up(i) right_patch_up_[(i)-1]
#define right_patch_down(i) right_patch_down_[(i)-1]

#define interaction_matrix(ipatch,jpatch) interaction_matrix_[indx2f(ipatch,jpatch,npatch)]
/*
 input:
 keep_left_states[] is length (target_up+1)*(target_down+1)
 keep_right_states[] is length (target_up+1)*(target_down+1)

 output:
 return npatch
 left_patch_size[] is length (target_up+1)*(target_down+1)
 right_patch_size[] is length (target_up+1)*(target_down+1)
*/

  assert( (bincoeff(10, 3 ) == 120) && (bincoeff(10,7) == bincoeff(10,3)) );

int nleft_up = 0;
int nleft_down = 0;
int nright_up = 0;
int nright_down = 0;

int npatch = 0;



/*
function [left_patch_size,right_patch_size] =  ...
          gen_patches_comb( left_size, right_size, ...
		       target_up, target_down, ...
                       keep_left_states, keep_right_states)
% ------------------------------------------------
%  [left_patch_size,right_patch_size] =  ...
%          gen_patches_comb( left_size, right_size, ...
%		       target_up, target_down, ...
%                      keep_left_states, keep_right_states)
%
%  estimate the size of patches based on 
%  combinatorial arguments
% ------------------------------------------------
*/

int max_left_up    = MIN( left_size,  target_up);
int max_left_down  = MIN( left_size,  target_down);
int max_right_up   = MIN( right_size, target_up);
int max_right_down = MIN( right_size, target_down);


/*
% ---------------------------------------------------------
% setup logical mask to identify valid (nleft_up,nleft_down)
% and valid matching (nright_up,nright_down)
% ---------------------------------------------------------
isvalid_left = zeros(1+max_left_up,1+max_left_down);
isvalid_right = zeros(1+max_right_up,1+max_right_down);
*/

int isvalid_left_[ (1+max_left_up)*(1+max_left_down)];
int isvalid_right_[ (1+max_right_up)*(1+max_right_down)];

#define isvalid_left(i,j)  isvalid_left_[ indx2f(i,j,(1+max_left_up)) ]
#define isvalid_right(i,j)  isvalid_right_[ indx2f(i,j,(1+max_right_up)) ]
  {
  for(nleft_down=0; nleft_down <= max_left_down; nleft_down++) {
  for(nleft_up=0;   nleft_up   <= max_left_up;   nleft_up++) {
       isvalid_left(1+nleft_up,1+nleft_down) = 0;
       };
       };

  for(nright_down=0; nright_down <= max_right_down; nright_down++) {
  for(nright_up=0;   nright_up   <= max_right_up;   nright_up++) {
       isvalid_right(1+nright_up,1+nright_down) = 0;
       };
       };
   }



/*
for nleft_up=0:max_left_up,
for nleft_down=0:max_left_down,
*/
for(nleft_down=0; nleft_down <= max_left_down; nleft_down++) {
for(nleft_up=0;   nleft_up   <= max_left_up;   nleft_up++) {
 
    int nright_up = target_up - nleft_up;
    int nright_down = target_down - nleft_down;
 
    int isvalid = (0 <= nright_up) && (nright_up <= max_right_up) && 
                  (0 <= nright_down) && (nright_down <= max_right_down);
    isvalid_left(1+nleft_up,1+nleft_down) = isvalid ? 1 : 0;
    isvalid_right(1+nright_up,1+nright_down) = isvalid ? 1 : 0;
 };
 };





/*
--------------------------------
npatch = sum( isvalid_left(:) );
--------------------------------
*/

  {
  int isum = 0;
  for(nleft_down=0; nleft_down <= max_left_down; nleft_down++) {
  for(nleft_up=0;   nleft_up   <= max_left_up;   nleft_up++) {
     int isvalid = isvalid_left( 1+nleft_up,1+nleft_down);
     if (isvalid) {  
        isum  = isum + 1;
        };
     };
     };
   npatch = isum;

   if (idebug >= 1) {
    printf("gen_patches_comb: isum=%d\n", isum );
    };

  }


/*
left_patch_size  = zeros(npatch,1);
right_patch_size = zeros(npatch,1);
*/





/*
left_ways  = zeros( max_left_up+1, max_left_down+1);
right_ways = zeros( max_right_up+1, max_right_down+1);
*/

double left_ways_[(max_left_up+1)*(max_left_down+1)];
double right_ways_[(max_right_up+1)*(max_right_down+1)];

#define left_ways(i,j) left_ways_[ indx2f(i,j,(1+max_left_up)) ]
#define right_ways(i,j) right_ways_[ indx2f(i,j, (1+max_right_up)) ]
  {
  for(nleft_down=0; nleft_down <= max_left_down; nleft_down++) {
  for(nleft_up=0; nleft_up <= max_left_up; nleft_up++) {
     left_ways(1+nleft_up,1+nleft_down) = 0;
     };
     };

  for(nright_down=0; nright_down <= max_right_down; nright_down++) {
  for(nright_up=0; nright_up <= max_right_up; nright_up++) {
     right_ways(1+nright_up,1+nright_down) = 0;
     };
     };
   }

  

/*
% ---------------------------------------------
% vectorized way to compute
% Cmat_left(1+nleft_up,1+nleft_down) = ...
%               nchoosek(left_size,nleft_up) *  ...
%               nchoosek(left_size,nleft_down)
% ---------------------------------------------
*/

/*
nleft_up     = 0:max_left_up;
nleft_down   = 0:max_left_down;
C_nleft_up   = bincoeff( left_size, nleft_up);
C_nleft_down = bincoeff( left_size, nleft_down);
*/

double C_nleft_up_[max_left_up+1];
double C_nleft_down_[max_left_down+1];

#define C_nleft_up(i) C_nleft_up_[(i)-1]
#define C_nleft_down(i) C_nleft_down_[(i)-1]

 {
 int nleft_up = 0;
 int nleft_down = 0;

 for(nleft_up = 0; nleft_up <= max_left_up; nleft_up++) {
   C_nleft_up(1+nleft_up) = bincoeff( left_size, nleft_up );
   };
 
 for(nleft_down = 0; nleft_down <= max_left_down; nleft_down++) {
    C_nleft_down(1+nleft_down) = bincoeff( left_size, nleft_down );
    };
 
  }


/*
%   ------------------------------------------------------
%   Cmat_left = C_nleft_up(:) * transpose(C_nleft_down(:));
%   ------------------------------------------------------
*/

double Cmat_left_[ (1+max_left_up)*(1+max_left_down)];
#define Cmat_left(i,j) Cmat_left_[ indx2f(i,j,(1+max_left_up)) ]

  {
  int nleft_up = 0;
  int nleft_down = 0;

  for( nleft_down = 0; nleft_down <= max_left_down; nleft_down++) {
  for( nleft_up = 0; nleft_up <= max_left_up; nleft_up++) {
     Cmat_left( 1+nleft_up,1+nleft_down) = C_nleft_up(1+nleft_up) * 
                                           C_nleft_down(1+nleft_down);
     };
     };
     
  }


/*
for nleft_down=0:max_left_down,
for nleft_up=0:max_left_up,
*/

  {
  int nleft_up = 0;
  int nleft_down = 0;

  for(nleft_down=0; nleft_down <= max_left_down; nleft_down++) {
  for(nleft_up=0; nleft_up <= max_left_up; nleft_up++) {

   if (isvalid_left(1+nleft_up,1+nleft_down)) {
/*
%    ---------------------------------
%    note, beware of possible overflow
%    or loss of accuracy
%    since nchoosek(n,k) can quickly become
%    very large numbers
%
%    for example, nchoosek(144,72) = 1.4802*10^42
%    ---------------------------------
*/
     left_ways(1+nleft_up,1+nleft_down) =  
         C_nleft_up(1+nleft_up) * 
         C_nleft_down(1+nleft_down);
     }
   };
   };
  }



/*
% ---------------------------------------------
% vectorized way to compute
% Cmat_right(1+nright_up,1+nright_down) = ...
%               nchoosek(right_size,nright_up) *  ...
%               nchoosek(right_size,nright_down)
% ---------------------------------------------
*/

/*
nright_up     = 0:max_right_up;
nright_down   = 0:max_right_down;
C_nright_up   = bincoeff( right_size, nright_up);
C_nright_down = bincoeff( right_size, nright_down);
*/

double C_nright_up_[ max_right_up+1];
double C_nright_down_[ max_right_down+1];

#define C_nright_up(i) C_nright_up_[(i)-1]
#define C_nright_down(i) C_nright_down_[(i)-1]

  {
  int nright_up = 0;
  int nright_down = 0;

  for(nright_up=0; nright_up <= max_right_up; nright_up++) {
     C_nright_up(1+nright_up) = bincoeff( right_size, nright_up );
     };

  for(nright_down=0; nright_down <= max_right_down; nright_down++) {
     C_nright_down(1+nright_down) = bincoeff( right_size, nright_down );
     };
  };




/*
%   ----------------------------------------------------------
%   Cmat_right = C_nright_up(:) * transpose(C_nright_down(:));
%   ----------------------------------------------------------
*/

double Cmat_right_[ (1+max_right_up)*(1+max_right_down)];
#define Cmat_right(i,j) Cmat_right_[ indx2f(i,j, (1+max_right_up)) ]


/*
for nright_down=0:max_right_down,
for nright_up = 0:max_right_up;
*/

  for(nright_down=0; nright_down <= max_right_down; nright_down++) {
  for(nright_up=0; nright_up <= max_right_up; nright_up++) {
     Cmat_right( 1+nright_up,1+nright_down) = C_nright_up(1+nright_up) * 
                                              C_nright_down(1+nright_down);
  };
  };
     



  {
  int nright_up = 0;
  int nright_down = 0;

  for(nright_down=0; nright_down <= max_right_down; nright_down++) {
  for(nright_up=0; nright_up <= max_right_up; nright_up++) {

   if (isvalid_right(1+nright_up,1+nright_down)) {
/*
%    ---------------------------------
%    note, beware of possible overflow
%    or loss of accuracy
%    since nchoosek(n,k) can quickly become
%    very large numbers
%
%    for example, nchoosek(144,72) = 1.4802*10^42
%    ---------------------------------
*/
     right_ways(1+nright_up,1+nright_down) =  
         C_nright_up(1+nright_up) * 
         C_nright_down(1+nright_down);
     }
   };
   };
  }

/*
total_left_ways = sum( left_ways(:) );
total_right_ways = sum( right_ways(:) );
*/

double total_left_ways = 0;
double total_right_ways = 0;
  {

  for(nleft_down=0; nleft_down <= max_left_down; nleft_down++) {
  for(nleft_up=0; nleft_up <= max_left_up; nleft_up++) {
    total_left_ways += left_ways(1+nleft_up,1+nleft_down);
    };
    };

  for(nright_down=0; nright_down <= max_right_down; nright_down++) {
  for(nright_up=0; nright_up <= max_right_up; nright_up++) {
    total_right_ways += right_ways(1+nright_up,1+nright_down);
    };
    };

   if (idebug >= 1) {
      printf("gen_patches_comb: total_left_ways=%lf, total_right_ways=%lf\n",
                                total_left_ways,     total_right_ways );
      };

   }

/*
% --------------------------------------------------------------------------------------
%
% set use_ceil = 1  to generate more patches
% set use_ceil = 0  to generate fewer patches
%
% --------------------------------------------------------------------------------------
% use ceiling() function  would generate at least 1 row per valid (nleft_up,nleft_down)
% this may generate more patches
%
% use round() function  would remove many unlikely or small states
% this may generate fewer patches
%
% using more patches would spread the rows across more patches and
% may perform less work overall
% --------------------------------------------------------------------------------------
*/

int use_ceil = 0;
/*
-------------------------------
if (use_ceil){
  left_ways = ceil((left_ways/total_left_ways) * keep_left_states );
  right_ways = ceil((right_ways/total_right_ways) * keep_right_states);

  }
else {
  left_ways = round((left_ways/total_left_ways) * keep_left_states );
  right_ways = round((right_ways/total_right_ways) * keep_right_states);
};
-------------------------------
*/

  for(nleft_up=0; nleft_up <= max_left_up; nleft_up++) {
  for(nleft_down=0; nleft_down <= max_left_down; nleft_down++) {
     double frac = ( left_ways(1+nleft_up,1+nleft_down)/total_left_ways * keep_left_states );

     left_ways( 1+nleft_up,1+nleft_down) =  (use_ceil) ? ceil( frac ) : round(frac);
     };
     };

  for(nright_up=0; nright_up <= max_right_up; nright_up++) {
  for(nright_down=0; nright_down <= max_right_down; nright_down++) {
     double frac = ( right_ways(1+nright_up,1+nright_down)/total_right_ways * keep_right_states );

     right_ways( 1+nright_up,1+nright_down) =  (use_ceil) ? ceil( frac ) : round(frac);
     };
     };



npatch = 0;

/*
for nleft_down=0:max_left_down,
for nleft_up=0:max_left_up,
*/

  for (nleft_down=0; nleft_down <= max_left_down; nleft_down++) {
  for (nleft_up=0; nleft_up <= max_left_up; nleft_up++) {
    nright_up = target_up - nleft_up;
    nright_down = target_down - nleft_down;
    
    
    if (isvalid_left(1+nleft_up,1+nleft_down)) {
     int lsize = left_ways(1+nleft_up,1+nleft_down);
     int rsize = right_ways(1+nright_up,1+nright_down);
  
     int has_work = (lsize >= 1) && (rsize >= 1);
     if (has_work) {
       npatch = npatch + 1;
       left_patch_size(npatch) = lsize;
       left_patch_up(npatch) = nleft_up;
       left_patch_down(npatch) = nleft_down;

       right_patch_size(npatch) = rsize;
       right_patch_up(npatch) = nright_up;
       right_patch_down(npatch) = nright_down;
     };
  
    };
  
  };
  };


         

/*
% -------------------------
% keep only non-zero patches
% -------------------------
*/

/*
idx_nonzero = find( (left_patch_size > 0) & (right_patch_size > 0) );
left_patch_size = left_patch_size(idx_nonzero);
right_patch_size = right_patch_size(idx_nonzero);
*/
  {
  int ipatch = 0;
  for(ipatch=1; ipatch <= npatch; ipatch++) {
    int isvalid = (left_patch_size(ipatch) > 0) &&
                  (right_patch_size(ipatch) > 0);
    assert( isvalid );
    };
  }



/*
    % ---------------------------------------------------------------
    % determine the sparse interaction matrix between (ipatch,jpatch)
    % ---------------------------------------------------------------
    interaction_matrix = zeros(npatch,npatch);
    for jpatch=1:npatch,
     nleft_up_j = patch2nleft_up(jpatch);
     nleft_down_j = patch2nleft_down(jpatch);
     for ipatch=1:npatch,
       nleft_up_i = patch2nleft_up(ipatch);
       nleft_down_i = patch2nleft_down(ipatch);
    
       njump = abs(nleft_up_i-nleft_up_j) + ...
               abs(nleft_down_i-nleft_down_j);
       interaction_matrix(ipatch,jpatch) = (njump <= 1);
     end;
    end;
*/
{
 int ipatch = 0;
 int jpatch = 0;
 
 for(jpatch=1; jpatch <= npatch; jpatch++) {
 for(ipatch=1; ipatch <= npatch; ipatch++) {
   int nleft_up_j = left_patch_up(jpatch);
   int nleft_down_j = left_patch_down(jpatch);

   int nleft_up_i = left_patch_up(ipatch);
   int nleft_down_i = left_patch_down(ipatch);

   int njump = ABS( nleft_up_i - nleft_up_j) + 
               ABS( nleft_down_i - nleft_down_j );
  
   interaction_matrix(ipatch,jpatch) = (njump <= 1);
   };
   };
}


  return( npatch );
}
    



