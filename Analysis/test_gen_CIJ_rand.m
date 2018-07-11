% ----------------------------------
% simple script to test gen_CIJ_rand
% ----------------------------------
addpath("../MiniApp");

total_sites = 32;
MaxKeep = 2000;

% ------------------
% assume at mid-point
% ------------------
nleft_sites = ceil( total_sites/2);
nright_sites = total_sites - nleft_sites;
ntarget_up = ceil( total_sites/2);
ntarget_down = ceil( total_sites/2);

% ---------------------------
% assume left side is growing
% ---------------------------
nleft_states = 4*MaxKeep;
nright_states = MaxKeep;

[left_patch_sizes,right_patch_sizes] = gen_patches_comb( nleft_sites, nright_sites, ...
                                                    ntarget_up, ntarget_down, ...
                                                    nleft_states, nright_states);
noperators = 4;
[veclen,nelem,gflops] = gen_CIJ_rand( left_patch_sizes, right_patch_sizes, noperators );

disp(sprintf('global size of system %g ', veclen ));
disp(sprintf('entries in CIJ %g ', nelem ));

% --------------------------
% test calling apply_Htarget
% --------------------------

X = rand(veclen,1);
tic();
Y = apply_Htarget( X );
elapsed_time = toc();
disp(sprintf('time for apply_Htarget is %g sec (%g gflops/sec) ', ...
              elapsed_time,   gflops/elapsed_time));
