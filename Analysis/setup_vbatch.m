function [Abatch,Bbatch,left_patch_start,right_patch_start,xy_patch_start] = ...
        setup_vbatch(noperator,npatches, left_patch_size, right_patch_size, ...
        interaction_matrix);
% [Abatch,Bbatch,left_patch_start,right_patch_start,xy_patch_start] = ...
%           setup_vbatch(noperator,npatches, left_patch_size, right_patch_size, ...
%           interaction_matrix);
%
% setup Abatch and Bbatch arrays and other data structures
%

ialign = 32;
left_max_state = sum(left_patch_size);
right_max_state = sum(right_patch_size);

nrow_Abatch = left_max_state;
ncol_Abatch = left_max_state * noperator;
ld_Abatch = ialign * ceil( left_max_state/ ialign );
size_Abatch = ld_Abatch  * ncol_Abatch;

nrow_Bbatch = right_max_state;
ncol_Bbatch = right_max_state * noperator;
ld_Bbatch = ialign * ceil( right_max_state/ialign );
size_Bbatch = ld_Bbatch * ncol_Bbatch;

Abatch =  zeros( ld_Abatch, ncol_Abatch );
Bbatch =  zeros( ld_Bbatch, ncol_Bbatch );

left_patch_start = cumsum( [1; left_patch_size(:)] );
right_patch_start = cumsum( [1; right_patch_size(:)] );

npatches = length(left_patch_size);
xy_patch_size = zeros(npatches,1);
for ipatch=1:npatches,
 nrowX = right_patch_size(ipatch);
 ncolX = left_patch_size(ipatch);
 xy_patch_size(ipatch) = nrowX * ncolX;
end;
xy_patch_start = cumsum( [1; xy_patch_size(:)]);

for i=1:nrow_Abatch,
for j=1:ncol_Abatch,
  ipatch = i;
  jpatch = j;
  has_work = interaction_matrix(ipatch,jpatch);
  if (has_work), 
    Abatch( i, j) = (i+(j-1)*nrow_Abatch) ./ (nrow_Abatch*ncol_Abatch);
  end;
end;
end;

for i=1:nrow_Bbatch,
for j=1:ncol_Bbatch,
  ipatch = i;
  jpatch = j;
  has_work = interaction_matrix(ipatch,jpatch);
  if (has_work),
    Bbatch(i,j) = (i+j) ./ (nrow_Bbatch * ncol_Bbatch);
  end;
end;
end;
