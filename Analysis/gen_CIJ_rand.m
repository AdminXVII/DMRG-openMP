function [veclen, nelem,gflops] = gen_CIJ_rand( left_patch_sizes, right_patch_sizes, noperators)
%
% [veclen, nelem,gflops] = gen_CIJ_rand( left_patch_sizes, right_patch_sizes, noperators)
%
% generate a CIJ data structure suitable for apply_Htarget
% but filled with random data in dense matrices
%
global CIJ;
global Lindex_patch;
global Rindex_patch;

CIJ = [];
Lindex_patch = [];
Rindex_patch = [];


% ---------------------------------------------------------
% option to add contribution of kron(HL,eye) + kron(eye,HR)
% ---------------------------------------------------------
add_identity = 0;

npatches = length(left_patch_sizes);
isok = (length(left_patch_sizes) == length(right_patch_sizes));
if (not(isok)),
  error(sprintf('gen_CIJ_rand: invalid number of patches (%d,%d) ', ...
                length(left_patch_sizes), length(right_patch_sizes) ) );
  return;
end;

% ----------------------------------------------
% setup Lindex_patch{:}, array of index vectors
% ----------------------------------------------
istart = 1;
for ipatch=1:npatches,
    iend = istart + left_patch_sizes(ipatch)-1;
    Lindex_patch{ipatch} = istart:iend;
    istart = iend + 1;
end;

% ----------------------------------------------
% setup Rindex_patch{:}, array of index vectors
% ----------------------------------------------
istart = 1;
for ipatch=1:npatches,
   iend = istart + right_patch_sizes(ipatch)-1;
   Rindex_patch{ipatch} = istart:iend;
   istart = iend + 1;
end;


% ---------
% setup CIJ
% ---------
for ipatch=1:npatches,
for jpatch=ipatch:npatches,
   
   left_nrow = left_patch_sizes(ipatch);
   left_ncol = left_patch_sizes(jpatch);

   right_nrow = right_patch_sizes(ipatch);
   right_ncol = right_patch_sizes(jpatch);

   Ak = [];
   Bk = [];
   for iop=1:noperators,
      Ak{iop} = rand( left_nrow, left_ncol );
      Bk{iop} = rand( right_nrow, right_ncol);
   end;

   
   CIJ{ipatch,jpatch}.Ak = Ak;
   CIJ{ipatch,jpatch}.Bk = Bk;

   is_diagonal = (ipatch == jpatch);
   if (not(is_diagonal)),
      % ----------------
      % enforce symmetry   
      % CIJ{ipatch,jpatch} = CIJ{jpatch,ipatch}'
      % ----------------

      Ak_trans = [];
      Bk_trans = [];

      for iop=1:length(Ak),
         Ak_trans{iop} = Ak{iop}';
         Bk_trans{iop} = Bk{iop}';
      end;

      CIJ{jpatch,ipatch}.Ak = Ak_trans;
      CIJ{jpatch,ipatch}.Bk = Bk_trans;
   else
    % -----------
    % is diagonal
    % -----------
    if (add_identity),
     % -----------------------------------------------
     % add contribution of kron(HL,eye) + kron(eye,HR)
     % note the apply_Htarget might add code to detect this
     % special case of identity matrix
     % -----------------------------------------------
     nops = length( CIJ{ipatch,jpatch}.Ak );

     % ----------------------------
     % contribution of kron(HL,eye)
     % ----------------------------
     CIJ{ipatch,jpatch}.Ak{nops+1} = eye( left_nrow, left_ncol);
     CIJ{ipatch,jpatch}.Bk{nops+1} = rand( right_nrow, right_ncol);

     % ----------------------------
     % contribution of kron(eye,HR)
     % ----------------------------
     CIJ{ipatch,jpatch}.Ak{nops+2} = rand( left_nrow, left_ncol);
     CIJ{ipatch,jpatch}.Bk{nops+2} = eye(  right_nrow, right_ncol);
    end;
   end;

end;
end;

% -------------------------------------------------------
% return vector length and total number of entries in CIJ
% -------------------------------------------------------
veclen = sum( left_patch_sizes(1:npatches) .* right_patch_sizes(1:npatches) );

nelem = 0;
gflops = 0;
for jpatch=1:npatches,
for ipatch=1:npatches,
   Ak = CIJ{ipatch,jpatch}.Ak;
   Bk = CIJ{ipatch,jpatch}.Bk;


   nop = length(Ak);
   nelem = nelem + prod(size(Ak{1})) * nop;
   nelem = nelem + prod(size(Bk{1})) * nop;

   [flops,imethod] = kron_total_flops( Ak{1}, Bk{1} );
   gflops = gflops + flops * nop;
end;
end;
  

gflops = gflops / 10^9;
