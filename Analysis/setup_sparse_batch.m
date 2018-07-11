function [state] = setup_sparse_batch( noperator, npatches, ...
                     left_patch_size, right_patch_size, ...
                     gnnz_A, gnnz_B)
% ----------------------------------------------------
% [state] = setup_sparse_batch( noperator, npatches, ...
%                     left_patch_size, right_patch_size, ...
%                     gnnz_A, gnnz_B)
% ----------------------------------------------------
 
idebug = 1;

% ---------------------
% setup data structures
% ---------------------
xy_patch_size = left_patch_size .* right_patch_size;
xy_patch_start = cumsum( [1; xy_patch_size(:)]);

nC = zeros(npatches,npatches);
for ioperator=1:noperator,
for jpatch=1:npatches,
for ipatch=1:npatches,
  is_zero_A = (gnnz_A(ipatch,jpatch,ioperator) <= 0);
  is_zero_B = (gnnz_B(ipatch,jpatch,ioperator) <= 0);
  if (is_zero_A || is_zero_A),
     % ---------------
     % no contribution
     % ---------------
  else
     nC(ipatch,jpatch) = nC(ipatch,jpatch) + 1;
  end;
end;
end;
end;

max_nC = max( nC(:) );

% ---------------------------------------
% construct compressed sparse row storage
% ---------------------------------------
ndeg = sum( (nC(1:npatches,1:npatches) >= 1), 2 );
xadj = cumsum( [1; ndeg(:)]);
jpatchC = zeros( sum(ndeg),1);
jnC = zeros( sum(ndeg),1);
ndeg = zeros(npatches,1);
for ipatch=1:npatches,
   for jpatch=1:npatches,
     if (nC(ipatch,jpatch) >= 1),
        ip = xadj(ipatch) + ndeg(ipatch);
        jnC(ip) = nC(ipatch,jpatch);
        jpatchC(ip) = jpatch;
        ndeg(ipatch) = ndeg(ipatch) + 1;
     end;
   end;
end;

% ------------
% double check
% ------------
ndeg = sum( nC(1:npatches,1:npatches) >= 1, 2 );
for ipatch=1:npatches,
   istart = xadj(ipatch);
   iend = xadj(ipatch+1)-1;
   isok = ((iend-istart+1) == ndeg(ipatch));
   if (~isok),
     disp(sprintf('setup_sparse_batch:ipatch=%d,istart=%d,iend=%d,ndeg(ipatch)=%d', ...
                   ipatch,   istart,   iend,   ndeg(ipatch) ));
   end;

   for k=istart:iend
     jpatch = jpatchC(k);
     isok = (jnC(k) == nC(ipatch,jpatch));
     if (~isok),
        disp(sprintf('setup_sparse_batch: ipatch=%d,jpatch=%d,k=%d,jnC(k)=%d,nC(ipatch,jpatch)=%d', ...
                      ipatch,   jpatch,   k,   jnC(k),   nC(ipatch,jpatch) ));
     end;
   end;
end;


  


% -----------------------------------------------
% setup Abatch, Bbatch as sparse data structures
% -----------------------------------------------

BXbatch_ncols = zeros(npatches,npatches);
Bbatch_ncols = zeros(npatches,npatches);
Abatch_ncols = zeros(npatches,npatches);

for ipatch=1:npatches,
for jpatch=1:npatches,

  nrowA = left_patch_size(ipatch);
  ncolA = left_patch_size(jpatch);

  nrowB = right_patch_size(ipatch);
  ncolB = right_patch_size(jpatch);
  
  nrowX = ncolB;
  ncolX = ncolA;

  nrowBX = nrowB;
  ncolBX = ncolX;
 
  nconnection = nC(ipatch,jpatch);

  Abatch_ncols(ipatch,jpatch)  = nconnection * ncolA;
  Bbatch_ncols(ipatch,jpatch)  = nconnection * ncolB;
  BXbatch_ncols(ipatch,jpatch) = nconnection * ncolBX;
end;
end;

% -----------------------------------
% setup starting locations of columns 
% on block row "ipatch"
% -----------------------------------
for ipatch=1:npatches,
   Abatch_start(ipatch,1:npatches) = ...
          cumsum( [1, Abatch_ncols(ipatch,1:(npatches-1))]); 
   Bbatch_start(ipatch,1:npatches) = ...
          cumsum( [1, Bbatch_ncols(ipatch,1:(npatches-1))]); 
   BXbatch_start(ipatch,1:npatches) = ...
          cumsum( [1, BXbatch_ncols(ipatch,1:(npatches-1))]); 
end;   

nA_rows = left_patch_size(1:npatches);
nB_rows = right_patch_size(1:npatches);
nBX_rows = nB_rows(1:npatches);

nA_cols = sum( Abatch_ncols, 2 );
nB_cols = sum( Bbatch_ncols, 2 );
nBX_cols = sum( BXbatch_ncols, 2 );





if (idebug >= 1),
   total_gAbatch  = sum( nA_rows(1:npatches)  .* nA_cols(1:npatches) );
   total_gBbatch  = sum( nB_rows(1:npatches)  .* nB_cols(1:npatches) );
   total_gBXbatch = sum( nBX_rows(1:npatches) .* nBX_cols(1:npatches) );

   disp(sprintf('total_gAbatch=%g, total_gBbatch=%g, total_gBXbatch=%g', ...
                 total_gAbatch,    total_gBbatch,    total_gBXbatch ));
end;

left_patch_start = cumsum( [1; left_patch_size(:)]);
right_patch_start = cumsum( [1; right_patch_size(:)]);

left_sum_size = sum( left_patch_size(:));
right_sum_size  = sum( right_patch_size(:));

if (idebug >= 1),
  disp(sprintf('noperator=%d, npatches=%d ', noperator,npatches));
  disp(sprintf('left_sum_size=%g, right_sum_size=%g', ...
                left_sum_size,    right_sum_size ));
end;

% -------------------------------
% setup values for Abatch, Bbatch
% -------------------------------
for ipatch=1:npatches,
  jpA = 1;
  jpB = 1;
  jpBX = 1;

  nrowA = left_patch_size(ipatch);
  nrowB = right_patch_size(ipatch);
  nrowBX = nrowB;

  sizeA = 0;
  sizeB = 0;
  sizeBX = 0;

  % ------------------------------------------------------------------
  % store Abatch as [A1 | A2 | ... ] (short and wide)
  % store Bbatch as [B1 | B2 | ... ] (short and wide)
  % store BXbatch as [BX1 | BX2 | ... ] (like Bbatch)
  % ------------------------------------------------------------------

  for jpatch=1:npatches,
    nconnection = nC(ipatch,jpatch);
    has_work = (nconnection >= 1);
    if (has_work),

     ncolA = left_patch_size(jpatch);
     ncolB = right_patch_size(jpatch);

     nrowX = ncolB;
     ncolX = ncolA;

     nrowBX = nrowB;
     ncolBX = ncolX;

     sizeA = sizeA + nconnection * ncolA;

     sizeB = sizeB + nconnection * ncolB;
     sizeBX = sizeBX + nconnection * ncolBX;
    end;
  end; % for jpatch

  % ------------
  % double check 
  % ------------
  isok = (sizeA == nA_cols(ipatch));
  if (~isok),
    error(sprintf('ipatch=%d, sizeA (%d) ~= nA_cols (%d) ', ...
                   ipatch,    sizeA,        nA_cols(ipatch) ));
    return;
  end;

  isok = (sizeB == nB_cols(ipatch));
  if (~isok),
    error(sprintf('ipatch=%d, sizeB (%d) ~= nB_cols (%d) ', ...
                   ipatch,    sizeB,        nB_cols(ipatch) ));
    return;
  end;

  isok = (sizeBX == nBX_cols(ipatch));
  if (~isok),
    error(sprintf('ipatch=%d, sizeBX (%d) ~= nBX_cols (%d) ', ...
                   ipatch,    sizeBX,        nBX_cols(ipatch) ));
    return;
  end;





  gAbatch{ipatch}.Abatch = zeros( nrowA, sizeA );
  gBbatch{ipatch}.Bbatch = zeros( nrowB, sizeB );
  gBXbatch{ipatch}.BXbatch  = zeros(nrowBX, sizeBX );

  for jpatch=1:npatches,

   nconnection = nC(ipatch,jpatch);
   has_work = (nconnection >= 1);
   if (has_work),
     ncolA = left_patch_size(jpatch);
     ncolB = right_patch_size(jpatch);

     nrowX = ncolB;
     ncolX = ncolA;

     nrowBX = ncolX;
     nrowBX = nrowB;

     Abatch  = zeros(nrowA,  nconnection * ncolA);
     Bbatch  = zeros(nrowB,  nconnection * ncolB);
     BXbatch = zeros(nrowBX, nconnection * ncolBX);

     for iconnection=1:nconnection,
         for ja=1:ncolA,
            ia=1:nrowA;
            gia = left_patch_start(ipatch) + (ia-1);
            gja = left_patch_start(jpatch) + (ja-1);

            aval = gia + (gja-1)*left_sum_size;
            aval = aval + (iconnection-1)*left_sum_size*left_sum_size;
            aval = aval / (max_nC * left_sum_size * left_sum_size);
            Amat(ia,ja) = aval;
         end;

         for jb=1:ncolB,
             ib=1:nrowB;
             gib = right_patch_start(ipatch) + (ib-1);
             gjb = right_patch_start(jpatch) + (jb-1);

             bval = gib + (gjb-1)*right_sum_size;
             bval = bval + (iconnection-1)*right_sum_size*right_sum_size;
             bval = bval / (max_nC * right_sum_size * right_sum_size );
             Bmat(ib,jb) = bval;
         end;

         i1 = 1;
         i2 = nrowB;

         j1 = 1 + (iconnection-1)*ncolB;
         j2 = j1 + ncolB-1;

         if (idebug >= 2),
           disp(sprintf('setup_sparse_batch: ipatch=%d,jpatch=%d,iconnection=%d,norm(Bmat)=%g', ...
                                             ipatch,   jpatch,   iconnection,   norm(Bmat,'fro')));
         end;

         Bbatch( i1:i2, j1:j2) = Bmat(1:nrowB,1:ncolB);


         i1 = 1;
         i2 = nrowA;

         j1 = 1 + (iconnection-1)*ncolA;
         j2 = j1 + ncolA-1;

         Abatch( i1:i2, j1:j2) = Amat(1:nrowA,1:ncolA);
      end; % for iconnection



    i1 = 1;
    i2 = nrowA;

    j1 = jpA;
    j2 = j1 + (ncolA*nconnection)-1;
    jpA = j2+1;

    ja_start = Abatch_start(ipatch,jpatch);
    isok = (j1 == ja_start);
    if (~isok),
     error(sprintf('ipatch=%d,jpatch=%d,j1=%d,ja_start=%d', ...
                    ipatch,   jpatch,   j1,   ja_start ));
     return;
    end;

    gAbatch{ipatch}.Abatch(i1:i2, j1:j2) = Abatch(1:nrowA, 1:(ncolA*nconnection));



    i1 = 1;
    i2 = nrowB;

    j1 = jpB;
    j2 = j1 + (ncolB*nconnection)-1;
    jpB = j2 + 1;

    jb_start = Bbatch_start(ipatch,jpatch);
    isok = (j1 == jb_start);
    if (~isok),
      error(sprintf('ipatch=%d,jpatch=%d,j1=%d, jb_start=%d', ...
                     ipatch,   jpatch,   j1,    jb_start ));
      return;
    end;

    gBbatch{ipatch}.Bbatch(i1:i2, j1:j2) = Bbatch(1:nrowB, 1:(ncolB*nconnection));

    if (idebug >= 2),
      disp(sprintf('ipatch=%d,jpatch=%d,j1=%d,j2=%d,norm(gBbatch.Batch(:,j1:j2))=%g', ...
                    ipatch,   jpatch,   j1,   j2, ...
                    norm(gBbatch{ipatch}.Bbatch(i1:i2,j1:j2),'fro') ));
    end;


    end; % end if
    
  end; % for jpatch


 end; % for ipatch

 if (idebug >= 2),
   for ipatch=1:npatches,
    disp(sprintf('ipatch=%d, norm(gBbatch{ipatch}.Bbatch)=%g ', ...
                  ipatch,    norm(gBbatch{ipatch}.Bbatch,'fro') ));
   end;
 end;
     

state.npatches = npatches;

state.left_patch_size = left_patch_size;
state.right_patch_size = right_patch_size;
state.left_patch_start = left_patch_start;
state.right_patch_start = right_patch_start;
    

state.xy_patch_size = xy_patch_size;
state.xy_patch_start = xy_patch_start;

   
state.nC = nC;

state.Abatch_ncols = Abatch_ncols;
state.Bbatch_ncols = Bbatch_ncols;
state.BXbatch_ncols = BXbatch_ncols;

state.Abatch_start = Abatch_start;
state.Bbatch_start = Bbatch_start;
state.BXbatch_start = BXbatch_start;

state.nA_rows = nA_rows;
state.nA_cols = nA_cols;

state.nB_rows = nB_rows;
state.nB_cols = nB_cols;

state.nBX_rows = nBX_rows;
state.nBX_cols = nBX_cols;

state.gAbatch = gAbatch;
state.gBbatch = gBbatch;
state.gBXbatch = gBXbatch;


% -----------------------------
% compressed sparse row storage
% -----------------------------
state.xadj = xadj;
state.jnC = jnC;
state.jpatchC = jpatchC;


