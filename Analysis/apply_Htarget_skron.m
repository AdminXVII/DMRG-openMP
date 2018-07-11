function yout = apply_Htarget_skron( xin, state )
% yout = apply_Htarget_skron( xin, state )
%
idebug = 1;

yout = zeros(size(xin));

npatches = state.npatches;
left_patch_size = state.left_patch_size;
right_patch_size = state.right_patch_size;

nC = state.nC;

if (idebug >= 1),
  Abatch_start = state.Abatch_start;
  Bbatch_start = state.Bbatch_start;
  BXbatch_start = state.BXbatch_start;
end;

xy_patch_size = left_patch_size(1:npatches) .* right_patch_size(1:npatches);
xy_patch_start = cumsum( [1;  xy_patch_size(:)]);

gAbatch = state.gAbatch;
gBbatch = state.gBbatch;


for ipatch=1:npatches,
  nrowA = left_patch_size(ipatch);
  nrowB = right_patch_size(ipatch);

  nrowY = nrowB;
  ncolY = nrowA;

  YI = zeros(nrowY,ncolY);

  Abatch = gAbatch{ipatch}.Abatch;
  Bbatch = gBbatch{ipatch}.Bbatch;

  ja_offset = 1;
  jb_offset = 1;
  for jpatch=1:npatches,
    nconnection = nC(ipatch,jpatch);
    has_work = (nconnection >= 1);
    if (has_work),

       ncolA = left_patch_size(jpatch);
       ncolB = right_patch_size(jpatch);

       nrowX = ncolB;
       ncolX = ncolA;

       ix1 = xy_patch_start(jpatch);
       ix2 = ix1 + (nrowX * ncolX)-1;
       XJ = reshape( xin( ix1:ix2 ), nrowX, ncolX );

       isok = (ja_offset == Abatch_start(ipatch,jpatch));
       if (~isok),
          disp(sprintf('apply_Htarget_skron: ipatch=%d,jpatch=%d,ja_offset=%d,Abatch_start=%d', ...
                        ipatch, jpatch, ja_offset, Abatch_start(ipatch,jpatch)));
       end;

       isok = (jb_offset == Bbatch_start(ipatch,jpatch));
       if (~isok),
          disp(sprintf('apply_Htarget_skron: ipatch=%d,jpatch=%d,jb_offset=%d,Bbatch_start=%d', ...
                        ipatch, jpatch, jb_offset, Bbatch_start(ipatch,jpatch)));
       end;
    

       ja1 = ja_offset;
       jb1 = jb_offset;
       for iconnection=1:nconnection,
           ja2 = ja1 + ncolA-1;
           jb2 = jb1 + ncolB-1;

           Amat = Abatch(1:nrowA, ja1:ja2);
           Bmat = Bbatch(1:nrowB, jb1:jb2);

           % -------------------------------------
           % Compute  YI = YI + kron(Amat,Bmat)*XJ
           % -------------------------------------
           YI(1:nrowY,1:ncolY) = YI(1:nrowY,1:ncolY)  + ...
              (Bmat(1:nrowB,1:ncolB) * XJ(1:nrowX,1:ncolX)) * ...
                 transpose(Amat(1:nrowA,1:ncolA));

           ja1 = ja2 + 1;
           jb1 = jb2 + 1;
       end;

       ja_offset = ja_offset + nconnection * ncolA;
       jb_offset = jb_offset + nconnection * ncolB;
     end; % if 
  end; 

  iy1 = xy_patch_start(ipatch);
  iy2 = iy1 + (nrowY * ncolY)-1;
  yout(iy1:iy2) = reshape( YI, (nrowY*ncolY),1);

end;


  
