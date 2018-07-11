function yout = apply_Htarget_sparse( xin, state )
% -----------------------------------------
% yout = apply_Htarget_sparse( xin, state )
% 
% perform matrix-vector multiply
%  yout = H * xin
% -----------------------------------------
idebug = 0;
%  --------------------
%  compute   BX = B * X
%  --------------------
npatches = state.npatches;
nC = state.nC;
left_patch_size = state.left_patch_size;
right_patch_size = state.right_patch_size;
gBbatch = state.gBbatch;
gAbatch = state.gAbatch;

if (idebug >= 2),
  disp(sprintf('apply_Htarget_sparse: '));
  for ipatch=1:npatches,
    disp(sprintf('ipatch=%d,norm(gBbatch{ipatch}.Bbatch)=%g', ...
                  ipatch, norm(gBbatch{ipatch}.Bbatch,'fro') ));
  end;
end;

yout = zeros(size(xin));
if (idebug >= 1),
  disp(sprintf('apply_Htarget_sparse: norm(xin,2)=%g', norm(xin,2) ));
end;

xy_patch_size = state.xy_patch_size;
xy_patch_start = state.xy_patch_start;


Bbatch_start = state.Bbatch_start;
BXbatch_start = state.BXbatch_start;
Abatch_start = state.Abatch_start;

nA_cols = state.nA_cols;
nBX_cols = state.nBX_cols;

for ipatch=1:npatches,
  for jpatch=1:npatches,
     nconnection = nC(ipatch,jpatch);
     has_work = (nconnection >= 1);
     if (has_work),
         nrowA = left_patch_size(ipatch);
         ncolA = left_patch_size(jpatch);

         nrowB = right_patch_size(ipatch);
         ncolB = right_patch_size(jpatch);

         nrowX = ncolB;
         ncolX = ncolA;

         nrowBX = nrowB;
         ncolBX = ncolX;

         ix1 = xy_patch_start(jpatch);
         ix2 = ix1 + xy_patch_size(jpatch) - 1;
         XJ = reshape( xin( ix1:ix2 ), nrowX,ncolX );


         % -----------------------------------------------------
         % compute   [B1*XJ | B2*XJ | ...] = [B1 | B2| ...] * XJ
         % -----------------------------------------------------

         bx_offset = BXbatch_start(ipatch,jpatch)-1;
         b_offset = Bbatch_start(ipatch,jpatch)-1;

         for iconnection=1:nconnection,
           jbx1 = bx_offset + 1 + (iconnection-1)*ncolBX;
           jbx2 = jbx1 + ncolBX-1;

           jb1 = b_offset + 1 + (iconnection-1)*ncolB;
           jb2 = jb1 + ncolB-1;
           state.gBXbatch{ipatch}.BXbatch(1:nrowB, jbx1:jbx2) = ...
              gBbatch{ipatch}.Bbatch(1:nrowB,jb1:jb2) * XJ;

         end;

     end;
   end; % for jpatch
   if (idebug >= 2),
     Bbatch = gBbatch{ipatch}.Bbatch;
     disp(sprintf('ipatch=%d, norm(Bbatch)=%g, size(Bbatch)=(%d,%d)', ...
                   ipatch,    norm(Bbatch,2), size(Bbatch,1), size(Bbatch,2)));
   end;
end; % for ipatch

% --------------------------
% perform YI = BX * trans(A)
% --------------------------

for ipatch=1:npatches,
   
   nrowA = left_patch_size(ipatch);
   nrowB = right_patch_size(ipatch);
   nrowBX = nrowB;

   nrowY = nrowB;
   ncolY = nrowA;

   isize = nA_cols(ipatch);
   isok = (isize == nBX_cols(ipatch));
   if (~isok),
     error(sprintf('ipatch=%d,nA_cols=%d,nBX_cols=%d', ...
                   ipatch, nA_cols(ipatch), nBX_cols(ipatch) ));
     return;
   end;

   BXbatch = state.gBXbatch{ipatch}.BXbatch;
   Abatch = state.gAbatch{ipatch}.Abatch;
   if (idebug >= 1 ),
    
     disp(sprintf('apply_Htarget_sparse: ipatch=%d,size(BXbatch)=(%d,%d), size(Abatch)=(%d,%d)', ...
                                         ipatch, ...
                                         size(BXbatch,1),size(BXbatch,2),  ...
                                         size(Abatch,1), size(Abatch,2) ));

     disp(sprintf('nrowA=%d,nrowB=%d,nA_cols=%d', ...
                   nrowA,nrowB, ...
                   nA_cols(ipatch)));

     disp(sprintf('norm(BXbatch)=%g, norm(Abatch)=%g', ...
                   norm(BXbatch,'fro'),  norm(Abatch,'fro') ));
   end;

   % --------------------------
   % form  YI = (BX)*transpose(A)
   % --------------------------
   YI = BXbatch(1:nrowBX,1:isize) * transpose( Abatch(1:nrowA,1:isize) );

   if (idebug >= 1),
     disp(sprintf('ipatch=%d,norm(YI,2)=%g, nrowY=%d,ncolY=%d', ...
                   ipatch,   norm(YI,2),    nrowY,   ncolY ));
   end;

   iy1 = xy_patch_start(ipatch);
   iy2 = iy1 + xy_patch_size(ipatch)-1;
   iylen = iy2-iy1+1;
   isok = (iylen == (nrowY*ncolY));
   if (~isok),
     error(sprintf('nrowY=%d,ncolY=%d,iy1=%d,iy2=%d,iylen=%d', ...
                   nrowY,   ncolY,   iy1,   iy2,   iylen ));
     return;
   end;
   yout(iy1:iy2) = reshape( YI, iylen,1);

end;
