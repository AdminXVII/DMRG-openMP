
Abatch;
Bbatch;
npatches = size(gnnz_A,1);
noperator = size(gnnz_A,3);

t1 = tic;
state = setup_sparse_batch( noperator, npatches, ...
                            left_patch_size, right_patch_size,
                            gnnz_A, gnnz_B );
disp(sprintf('setup_sparse_batch took %g sec', toc(t1)  ));
nC = state.nC;
dummy = write_batch_info( npatches, left_patch_size, right_patch_size,  nC );

xy_sum_size = sum( left_patch_size .* right_patch_size );
xin = reshape( 1:xy_sum_size,   xy_sum_size,1);

t1 = tic;
yout = apply_Htarget_sparse( xin, state );
disp(sprintf('apply_Htarget_sparse took %g sec', toc(t1)));

yout2 = apply_Htarget_skron( xin, state);
abserr = norm(yout-yout2);
relerr = abserr/max( norm(yout), norm(yout2) );
disp(sprintf('abserr=%g, relerr=%g from apply_Htarget_skron', ...
              abserr,    relerr ));

disp(sprintf('norm(yout)=%g', norm(yout) ));

avg_yout = sum(yout)/prod(size(yout));
max_yout = max( yout );
min_yout = min( yout );

sd_yout = sqrt( sum( (yout - avg_yout).*(yout - avg_yout) ) );

disp(sprintf('avg(yout)=%g, max(yout)=%g, min(yout)=%g sd_yout=%g', ...
              avg_yout,    max_yout,    min_yout, sd_yout ));


