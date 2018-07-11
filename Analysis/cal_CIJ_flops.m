function [flops_total,flops_CIJ] = cal_CIJ_flops( npatches, ...
                                                  left_patch_size, ...
                                                  right_patch_size,  ...
                                                  interaction_matrix )
%
% [flops_total,flops_CIJ] = cal_CIJ_flops( npatches, ...
%                                          left_patch_size, ...
%                                          right_patch_size, ...
%                                          interaction_matrix )
%
% estimate the work for just  a single kron(Ak,Bk) part
%

flops_CIJ = zeros(npatches,npatches);

for jpatch=1:npatches,
for ipatch=1:npatches,
 has_work = interaction_matrix(ipatch,jpatch);
 if (has_work),
  % --------------------------------------------------------
  % In cell CIJ(ipatch,jpatch)
  % Ak is  left_patch_size(ipatch) by left_patch_size(jpatch)
  % Bk is  right_patch_size(ipatch) by right_patch_size(jpatch)
  % --------------------------------------------------------
  nrowA = left_patch_size(ipatch);
  ncolA = left_patch_size(jpatch);

  nrowB = right_patch_size(ipatch);
  ncolB = right_patch_size(jpatch);

  flops_CIJ(ipatch,jpatch) = cal_kron_flops( nrowA,ncolA, nrowB,ncolB);
 end;
end;
end;

flops_total = sum( sum( flops_CIJ ) );
  


