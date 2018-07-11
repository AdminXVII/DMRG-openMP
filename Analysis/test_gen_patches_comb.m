% ---------------------------------------
% simple script to test gen_patches_comb
% ---------------------------------------
 left_size = 8;
 right_size = 10;
 target_up = (left_size + right_size)/2;
 target_down = target_up;
 keep_left_states = 5000;
 keep_right_states = 4*keep_left_states;
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
for ipatch=1:npatches,
   printf("ipatch=%d, left_patch_size=%d, right_patch_size=%d\n",
           ipatch, left_patch_size(ipatch), right_patch_size(ipatch) );
end;

nnz_interaction_matrix = nnz( interaction_matrix );
disp(sprintf('number of non-zeros in interaction_matrix is %g', ...
     nnz_interaction_matrix ));

figure(1);
clf;
spy( interaction_matrix );
title(sprintf('interaction matrix, npatches=%d',npatches));


