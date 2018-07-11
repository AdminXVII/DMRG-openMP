function retval = write_CIJ_shape( npatches, left_patch_size, right_patch_size, fname , vectorLength)
%
% retval = cal_CIJ_flops( npatches, left_patch_size, right_patch_size, fname )
%
% Write out the shape and dimensions of Ak and Bk entries in CIJ
%

fd1 = fopen(fname, 'w');
fprintf(fd1, "%d  %d\n", npatches, vectorLength);
for ipatch=1:npatches,
 for jpatch=1:npatches,
  % --------------------------------------------------------
  % In cell CIJ(ipatch,jpatch)
  % Ak is  left_patch_size(ipatch) by left_patch_size(jpatch)
  % Bk is  right_patch_size(ipatch) by right_patch_size(jpatch)
  % --------------------------------------------------------
  nrowA = left_patch_size(ipatch);
  ncolA = left_patch_size(jpatch);
  numbytes = fprintf(fd1, "%d %d A %d %d", ipatch, jpatch, nrowA, ncolA);
  nrowB = right_patch_size(ipatch);
  ncolB = right_patch_size(jpatch);
  numbytes = fprintf(fd1, " B %d %d\n", nrowB, ncolB);
 
end;
end;
fclose(fd1);
retval = 0;  


