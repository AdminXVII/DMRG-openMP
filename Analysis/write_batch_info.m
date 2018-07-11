function dummy = write_batch_info( npatches, left_patch_size, right_patch_size,  nC )
%
% dummy = write_batch_info( npatches, left_patch_size, right_patch_size,  nC )
%
dummy = 0;

fid = fopen('left_patch_size.txt','w');
fprintf(fid,'%d\n', npatches);
for ipatch=1:npatches,
  fprintf(fid,'%d\n', left_patch_size(ipatch));
end;
fclose(fid);

fid = fopen('right_patch_size.txt','w');
fprintf(fid,'%d\n', npatches);
for ipatch=1:npatches,
  fprintf(fid,'%d\n', right_patch_size(ipatch));
end;
fclose(fid);

fid = fopen('nC.txt','w');
fprintf(fid,'%d\n', npatches);
for ipatch=1:npatches,
for jpatch=1:npatches,
  if (nC(ipatch,jpatch) ~= 0),
    fprintf(fid,'%d %d %d\n', ipatch,jpatch,nC(ipatch,jpatch));
  end;
end;
end;
fclose(fid);

