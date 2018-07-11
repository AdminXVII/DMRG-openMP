% ---------------------------
% read in sizes and non-zeros
% ---------------------------
Abatch;
Bbatch;
npatches = size(gnnz_A,1);
noperator = size(gnnz_A,3);


for ioperator=1:noperator,
  nnz_A(ioperator) = nnz( gnnz_A(:,:,ioperator));
  nnz_B(ioperator) = nnz( gnnz_B(:,:,ioperator));
end;

do_plot = 0;
if (do_plot),
  plot( 1:noperator, nnz_A,'r', 1:noperator, nnz_B, 'b');
  legend('nnz_A', 'nnz_B');
  title(sprintf('npatches=%d, number of nonzero patches',npatches));
end;


max_left_size = sum( left_patch_size );
max_right_size = sum( right_patch_size );

disp(sprintf('npatches=%d, noperator=%d', npatches,noperator));
disp(sprintf('max_left_size=%d, max_right_size=%d', max_left_size, max_right_size));


pattern = (gnnz_A ~= 0) & (gnnz_B ~= 0);
all_pattern = pattern(:,:,1);
for ioperator=2:noperator,
   all_pattern = all_pattern | pattern(:,:,ioperator);
end;

disp(sprintf('nnz(all_pattern)=%d ', nnz(all_pattern)));

% --------------------------------------
% number of operators per (ipatch,jpatch)
% --------------------------------------
nC = zeros( npatches, npatches );
for ioperator=1:noperator,
for jpatch=1:npatches,
for ipatch=1:npatches,
    non_zero_A = (gnnz_A(ipatch,jpatch,ioperator) > 0);
    non_zero_B = (gnnz_B(ipatch,jpatch,ioperator) > 0);
    if (non_zero_A && non_zero_B),
       nC(ipatch,jpatch) = nC(ipatch,jpatch) + 1;
    end;
end;
end;
end;

disp('nC')
nC




% ----------------------
% compute sparsity ratio
% ----------------------
ratio_A = zeros([npatches,npatches,noperator]);
ratio_B = zeros([npatches,npatches,noperator]);
for ioperator=1:noperator,
 for jpatch=1:npatches,
  ipatch=1:npatches;
     nrowA = left_patch_size(ipatch);
     ncolA = left_patch_size(jpatch);
     ratio_A(ipatch,jpatch,ioperator) = gnnz_A(ipatch,jpatch,ioperator) ./ (nrowA * ncolA);

     nrowB = right_patch_size(ipatch);
     ncolB = right_patch_size(jpatch);
     ratio_B(ipatch,jpatch,ioperator) = gnnz_B(ipatch,jpatch,ioperator) ./ (nrowB * ncolB);
  end;
end;

ioperator = 1;
max_ratio_A = reshape( ratio_A(:,:,ioperator), [npatches,npatches]);
max_ratio_B = reshape( ratio_B(:,:,ioperator), [npatches,npatches]);
avg_ratio_A = reshape( ratio_A(:,:,ioperator), [npatches,npatches]);
avg_ratio_B = reshape( ratio_B(:,:,ioperator), [npatches,npatches]);
for ioperator=2:noperator,
  max_ratio_A = max( max_ratio_A, ratio_A(:,:,ioperator) );
  max_ratio_B = max( max_ratio_B, ratio_B(:,:,ioperator) );
  avg_ratio_A = avg_ratio_A + ratio_A(:,:,ioperator);
  avg_ratio_B = avg_ratio_B + ratio_B(:,:,ioperator);
end;
avg_ratio_A = avg_ratio_A / noperator;
avg_ratio_B = avg_ratio_B / noperator;


% -------------
% estimate work 
% -------------
flops_method_1 = zeros([npatches,npatches,noperator]);
flops_method_2 = zeros([npatches,npatches,noperator]);
flops_method_min = zeros([npatches,npatches,noperator]);

for ioperator=1:noperator,
  for jpatch=1:npatches,
    for ipatch=1:npatches,
       nrowA = left_patch_size(ipatch);
       ncolA = left_patch_size(jpatch);
       nrowB = right_patch_size(ipatch);
       ncolB = right_patch_size(jpatch);

       has_work = (gnnz_A(ipatch,jpatch,ioperator) > 0) & ...
                  (gnnz_B(ipatch,jpatch,ioperator) > 0);
       if (has_work),
         [flops_total, flops_method1, flops_method2] = ...
             cal_kron_flops( nrowA,nrowB, ncolA, ncolB );

         flops_min = min(flops_method1,flops_method2);
         flops_method_1(ipatch,jpatch,ioperator) = flops_method1;
         flops_method_2(ipatch,jpatch,ioperator) = flops_method2;
     
         flops_method_min(ipatch,jpatch,ioperator) = flops_min;
      end;
    end;
 end;
end;

total_flops_method_1 = sum( flops_method_1(:) );
total_flops_method_2 = sum( flops_method_2(:) );
total_flops_method_min = sum( flops_method_min(:) );
ratio1 = total_flops_method_1 / total_flops_method_min;
ratio2 = total_flops_method_2 / total_flops_method_min;

disp(sprintf('flops1=%g (%g), flops2=%g (%g), flops_min=%g ', ...
      total_flops_method_1, ratio1, ...
      total_flops_method_2, ratio2, ...
      total_flops_method_min ));


% ----------------------
% estimate sparse flops
% ----------------------

sflops_method_1 = zeros([npatches,npatches,noperator]);
sflops_method_2 = zeros([npatches,npatches,noperator]);
sflops_method_3 = zeros([npatches,npatches,noperator]);
sflops_method_min = zeros([npatches,npatches,noperator]);

for ioperator=1:noperator,
  for jpatch=1:npatches,
    for ipatch=1:npatches,
       nrowA = left_patch_size(ipatch);
       ncolA = left_patch_size(jpatch);
       nrowB = right_patch_size(ipatch);
       ncolB = right_patch_size(jpatch);
       nnzA = gnnz_A(ipatch,jpatch,noperator);
       nnzB = gnnz_B(ipatch,jpatch,noperator);


       has_work = (gnnz_A(ipatch,jpatch,ioperator) > 0) & ...
                  (gnnz_B(ipatch,jpatch,ioperator) > 0);
       if (has_work),
         [sflops_total, sflops_method1, sflops_method2, sflops_method3] = ...
             cal_kron_sflops( nrowA,nrowB, ncolA, ncolB, nnzA, nnzB );

         sflops_min = min( min(sflops_method1,sflops_method2), sflops_method3);

         sflops_method_1(ipatch,jpatch,ioperator) = sflops_method1;
         sflops_method_2(ipatch,jpatch,ioperator) = sflops_method2;
         sflops_method_3(ipatch,jpatch,ioperator) = sflops_method3;
     
         sflops_method_min(ipatch,jpatch,ioperator) = sflops_min;
      end;
    end;
 end;
end;

total_sflops_method_1 = sum( sflops_method_1(:) );
total_sflops_method_2 = sum( sflops_method_2(:) );
total_sflops_method_3 = sum( sflops_method_3(:) );
total_sflops_method_min = sum( sflops_method_min(:) );

ratio1 = total_sflops_method_1 / total_sflops_method_min;
ratio2 = total_sflops_method_2 / total_sflops_method_min;
ratio3 = total_sflops_method_3 / total_sflops_method_min;

disp(sprintf('sflops1=%g (%g), sflops2=%g (%g), sflops3=%g (%g), sflops_min=%g ', ...
      total_sflops_method_1, ratio1, ...
      total_sflops_method_2, ratio2, ...
      total_sflops_method_3, ratio3, ...
      total_sflops_method_min ));



