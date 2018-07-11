% ---------------------------------------
% simple script to test work in method 1 vs method 2
% ---------------------------------------
 left_size = 8;
 right_size = 10;
 target_up = (left_size + right_size)/2;
 target_down = target_up;
 max_keep_states = 5000;
 keep_left_states = max_keep_states;
 keep_right_states = 4*keep_left_states;
 max_patches = (1 + target_up) * (1 + target_down);


[left_patch_size,right_patch_size] =  ...
          gen_patches_comb( left_size, right_size, ...
		       target_up, target_down, ...
                       keep_left_states, keep_right_states);

npatches = length( left_patch_size );

use_sort = 0;
if (use_sort),
  left_patch_size = sort(left_patch_size);
  right_patch_size = sort(right_patch_size);
end;


flops_method1_all = zeros(npatches,npatches);
flops_method2_all = zeros(npatches,npatches);
space_method1_all = zeros(npatches,npatches);
space_method2_all = zeros(npatches,npatches);

for jpatch=1:npatches,
for ipatch=1:npatches,
  nrowA = left_patch_size(ipatch);
  ncolA = left_patch_size(jpatch);
  nrowB = right_patch_size(ipatch);
  ncolB = right_patch_size(jpatch);
  [flops_total,flops_method1, flops_method2] = cal_kron_flops(nrowA,nrowB,ncolA,ncolB);
 
  flops_method1_all(ipatch,jpatch) = flops_method1;
  flops_method2_all(ipatch,jpatch) = flops_method2;

  [space_method1, space_method2] = cal_kron_space(nrowA,nrowB,ncolA,ncolB);
  space_method1_all(ipatch,jpatch) = space_method1;
  space_method2_all(ipatch,jpatch) = space_method2;
  
end;
end;

size_X = sum( left_patch_size .* right_patch_size );
disp(sprintf('size of X vector = %g', size_X ));



min_total_flops = sum( min( flops_method1_all(:), flops_method2_all(:) ) );
max_total_flops = sum( max( flops_method1_all(:), flops_method2_all(:) ) );

disp(sprintf('min_total_flops=%g, max_total_flops=%g', ...
     min_total_flops, max_total_flops ));
disp(sprintf('total_flops_method1=%g, total_flops_method2=%g', ...
       sum( flops_method1_all(:) ),  sum( flops_method2_all(:) )  ));


min_total_space = sum( min( space_method1_all(:), space_method2_all(:)));
max_total_space = sum( max( space_method1_all(:), space_method2_all(:)));

disp(sprintf('min_total_space=%g, max_total_space=%g', ...
              min_total_space, max_total_space ));

disp(sprintf('total_space_method1=%g, total_space_method2=%g', ...
              sum( space_method1_all(:)),  sum( space_method2_all(:)) ));

use_method1 = (flops_method1_all <= flops_method2_all);
spy( use_method1 );


