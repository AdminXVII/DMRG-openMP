function [Y] = apply_Htarget_vbatch( noperator,npatches, ...
          left_patch_start,right_patch_start,xy_patch_start, ...
          Abatch, Bbatch, X, interaction_matrix );
%
%[Y] = apply_Htarget_vbatch( noperator,npatches, ...
%          left_patch_start,right_patch_start,xy_patch_start, ...
%          Abatch, Bbatch, X, interaction_matrix );
%
% perform kronecker product 
%
left_max_states = left_patch_start(npatches+1)-1;
right_max_states = right_patch_start(npatches+1)-1;

nrowA = left_max_states;
ncolA = nrowA;
nrowB = right_max_states;
ncolB = nrowB;

left_patch_size = left_patch_start(1 + (1:npatches)) - ...
                  left_patch_start( 1:npatches );
right_patch_size = right_patch_start(1 + (1:npatches)) - ...
                  right_patch_start( 1:npatches );

xy_size = xy_patch_start(npatches+1)-1;
Y = zeros( xy_size,1);
for ipatch=1:npatches,
  nrowY = right_patch_size(ipatch);
  ncolY = left_patch_size(ipatch);

  YI = zeros(nrowY,ncolY);
  for jpatch=1:npatches,
   has_work = interaction_matrix(ipatch,jpatch);
   if (has_work),
    j1 = xy_patch_start(jpatch);
    j2 = xy_patch_start(jpatch+1)-1;

    R1 = right_patch_start(jpatch);
    R2 = right_patch_start(jpatch+1)-1;
    L1 = left_patch_start(jpatch);
    L2 = left_patch_start(jpatch+1)-1;

    nrowX = R2-R1+1;
    ncolX = L2-L1+1;
    XJ =  reshape( X( j1:j2 ), nrowX, ncolX );

    ib1 = right_patch_start(ipatch);
    ib2 = right_patch_start(ipatch+1)-1;
    jb1 = right_patch_start(jpatch);
    jb2 = right_patch_start(jpatch+1)-1;

    ia1 = left_patch_start(ipatch);
    ia2 = left_patch_start(ipatch+1)-1;
    ja1 = left_patch_start(jpatch);
    ja2 = left_patch_start(jpatch+1)-1;

    for k=1:noperator,
       offsetA = (k-1) * ncolA;
       offsetB = (k-1) * ncolB;
       AIJ = Abatch( ia1:ia2,  offsetA + (ja1:ja2) );
       BIJ = Bbatch( ib1:ib2,  offsetB + (jb1:jb2) );
       YI = YI + (BIJ * XJ) * transpose( AIJ );
    end; % for k
   end; % if (has_work)
  end; % for jpatch

   i1 = xy_patch_start(ipatch);
   i2 = xy_patch_start(ipatch+1)-1;
   Y(i1:i2) = reshape( YI, (i2-i1+1),1);
end;
