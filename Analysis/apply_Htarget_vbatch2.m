function [Y] = apply_Htarget_vbatch2( noperator,npatches, ...
          left_patch_start,right_patch_start,xy_patch_start, ...
          Abatch, Bbatch, X );
%
%[Y] = apply_Htarget_vbatch2( noperator,npatches, ...
%          left_patch_start,right_patch_start,xy_patch_start, ...
%          Abatch, Bbatch, X );
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


  % 
  % ----------------------------------------------------------
  % form  BX(:, ja1:ja2) = B(:,jb1:jb2) * XJ(jb1:jb2, ja1:ja2)
  % ----------------------------------------------------------
  nrowBX = size(Bbatch,1); 
  ncolBX = size(Abatch,2);
  BX = zeros( nrowBX, ncolBX );

  for jpatch=1:npatches,
    j1 = xy_patch_start(jpatch);
    j2 = xy_patch_start(jpatch+1)-1;

    R1 = right_patch_start(jpatch);
    R2 = right_patch_start(jpatch+1)-1;
    L1 = left_patch_start(jpatch);
    L2 = left_patch_start(jpatch+1)-1;

    nrowX = R2-R1+1;
    ncolX = L2-L1+1;
    XJ =  reshape( X( j1:j2 ), nrowX, ncolX );


    for k=1:noperator,
       offsetB = (k-1) * ncolB;
       offsetBX = (k-1) * ncolA;

       BX(1:nrowBX, offsetBX + (L1:L2) ) = ...
             Bbatch(1:nrowBX, offsetB + (R1:R2)) * XJ(1:(R2-R1+1),1:(L2-L1+1));
    end;
   end;

   Y = zeros( xy_size,1);
   for ipatch=1:npatches,
     L1 = left_patch_start(ipatch);
     L2 = left_patch_start(ipatch+1)-1;
     R1 = right_patch_start(ipatch);
     R2 = right_patch_start(ipatch+1)-1;

     YI = BX(R1:R2, 1:ncolBX) * ...
              transpose(Abatch(L1:L2,1:ncolBX));
     i1 = xy_patch_start(ipatch);
     i2 = xy_patch_start(ipatch+1)-1;
     Y(i1:i2) = reshape(YI, (i2-i1+1),1);
    end;

   



