function [flops_sweep_total] = sim_sweep( lattice_size, max_keep_states, do_plot, do_write)
% [flops_sweep_total] = sim_sweep( lattice_size, max_keep_states, do_plot, do_write)
% 
% simple script to simulate the work for a sweep
%

do_sweep = 1;
%do_plot = 1;

target_up = round(lattice_size/2);
target_down = round(lattice_size/2);
flops_work = zeros( lattice_size-1,1);


flops_sweep_total = 0;

if (do_sweep),
  istart = 1;
  iend = (lattice_size-1);
else
  % ------------------
  % just for debugging
  % ------------------
  istart = round(lattice_size/2);
  iend = istart;
end;

vector_length = zeros(lattice_size,1);

isize=iend-istart+1;
for left_size = istart:iend,
    right_size = lattice_size - left_size;

    % ----------------------------
    % assume growing the left part
    % ----------------------------
    keep_left_states = min( 4^left_size, 4*max_keep_states );
    keep_right_states = min( 4^right_size, max_keep_states );

    [left_patch_size,right_patch_size, interaction_matrix] =  ...
          gen_patches_comb( left_size, right_size, ...
		       target_up, target_down, ...
                       keep_left_states, keep_right_states);



    npatches = length(left_patch_size);
    [flops_total,flops_CIJ] = cal_CIJ_flops( npatches, ...
                                             left_patch_size, ...
                                             right_patch_size, ...
                                             interaction_matrix );

    flops_sweep_total = flops_sweep_total + flops_total;

    flops_work(left_size) = flops_total;
    vector_length(left_size) = sum( left_patch_size(1:npatches) .* ...
                                    right_patch_size(1:npatches) );
    if (do_write),
      fname = sprintf("CIJ_shape_%d.txt", left_size); 
      retval = write_CIJ_shape( npatches, ...
                                      left_patch_size, ...
                                      right_patch_size, ...
                                      fname, ...
                                      vector_length(left_size) );
    end,
                                      
    if (left_size == round(lattice_size/2)),
       left_patch_size_save = left_patch_size;
       right_patch_size_save = right_patch_size;
    end;
end;

% --------------
% generate plots
% --------------

if (do_plot >= 1),
  figure(1);
  clf;

      left_patch_size = left_patch_size_save;
      right_patch_size = right_patch_size_save;
      npatches = length(left_patch_size);
      [flops_total,flops_CIJ] = cal_CIJ_flops(npatches,...
                                       left_patch_size,...
                                       right_patch_size);

      left_size = round(lattice_size/2);
      right_size = lattice_size  - left_size;

      subplot(2,2,1);
      mesh( flops_CIJ );
      title(sprintf('size=(%d,%d),npatches=%d,\n states=(%d,%d),flops total=%g', ...
                     left_size, right_size,   ...
                     npatches, ...
                     sum(left_patch_size_save),  ...
                     sum(right_patch_size_save), ...
                     flops_total ));
      

      subplot(2,2,2);

      use_histogram = 1;
      if (use_histogram),
       % -----------------------------------------------
       % adjust the value of bins based on actual values
       % -----------------------------------------------
       imax =  ceil( log10( max(flops_CIJ(:)) ) );
       imin = floor( log10( min(flops_CIJ(:)) ) );
       hist( log10(flops_CIJ(:) ), imin:imax );
       xlabel('work (log10)');
       ylabel('number of cells');
       title('histogram of work in cells at mid-point');
      else
        nn = prod(size(flops_CIJ));
        semilogy( 1:nn, sort(flops_CIJ(:))  );
        title('sorted profile of flops in CIJ cells');
      end;

      subplot(2,2,3);
      semilogy( 1:size(flops_CIJ,1), sum(flops_CIJ,2) );
      title('sum of work by row in CIJ');


  if (do_sweep),
     subplot(2,2,4);
     semilogy( istart:iend, flops_work(istart:iend) ,'.-');
     xlabel('left size');
     ylabel('flops');
     title(sprintf('total work over a sweep=%g,\n lattice size=%d',...
         flops_sweep_total, lattice_size));
  end;
print(sprintf('size%d_%d.png',lattice_size,max_keep_states),'-dpng');


figure(2);
clf;
semilogy( istart:iend, vector_length(istart:iend), 'b-d');
title('length of vector');
xlabel('left size');
print(sprintf('veclen%d_%d.png',lattice_size,max_keep_states),'-dpng');




end;

