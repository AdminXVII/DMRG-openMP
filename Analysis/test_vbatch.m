
% ---------------------------------------
% simple script to test test_vbatch
% ---------------------------------------
idebug = 1;
 noperator = 4;
 left_size = 8;
 right_size = 8;
 target_up = (left_size + right_size)/2;
 target_down = target_up;
 max_keep_states = 1024;

 % ---------------------------
 % assume left part is growing
 % ---------------------------
 keep_left_states = 4*max_keep_states;
 keep_right_states = max_keep_states;
 max_patches = (1 + target_up) * (1 + target_down);


[left_patch_size,right_patch_size,interaction_matrix] =  ...
          gen_patches_comb( left_size, right_size, ...
		       target_up, target_down, ...
                       keep_left_states, keep_right_states);

npatches = length( left_patch_size );


printf("left_size=%d, right_size=%d, target_up=%d, target_down=%d\n", 
         left_size,    right_size,    target_up,    target_down );
printf("keep_left_states=%d, keep_right_states=%d\n",
         keep_left_states,    keep_right_states );

printf("npatches=%d\n", npatches );
if (idebug >= 2),
  for ipatch=1:npatches,
   printf("ipatch=%d, left_patch_size=%d, right_patch_size=%d\n",
           ipatch, left_patch_size(ipatch), right_patch_size(ipatch) );
  end;
end;

[Abatch,Bbatch, left_patch_start, right_patch_start, xy_patch_start] = ...
         setup_vbatch(noperator,npatches, ...
                      left_patch_size, right_patch_size, ...
                      interaction_matrix);


xy_size = xy_patch_start(npatches+1)-1;
X = reshape( (1:xy_size)/xy_size, xy_size,1);
Y = zeros( xy_size,1);


time_Y = tic();
Y = apply_Htarget_vbatch( noperator,npatches, ...
        left_patch_start,right_patch_start,xy_patch_start, ...
        Abatch, Bbatch, X, interaction_matrix );                       
time_Y = toc(time_Y);

time_Y2 = tic();
Y2 = apply_Htarget_vbatch2( noperator,npatches, ...
        left_patch_start,right_patch_start,xy_patch_start, ...
        Abatch, Bbatch, X, interaction_matrix );                       
time_Y2 = toc(time_Y2);

disp(sprintf('time for Y=%g sec, time for Y2=%g sec', ...
                time_Y,          time_Y2 ));

relerr = norm(Y-Y2,1)/max( norm(Y,1), norm(Y2,1) );
disp(sprintf('relerr=%g, norm(Y-Y2)=%g, norm(Y)=%g, norm(Y2)=%g ', ...
              relerr,    norm(Y-Y2,1),  norm(Y,1),  norm(Y2,1) ));


Y_avg = sum(Y)/xy_size;
Y_sd = sqrt(sum( (Y - Y_avg) .* (Y - Y_avg) ));


disp(sprintf('Y_avg = %g, Y_sd = %g', ...
              Y_avg,      Y_sd ));





