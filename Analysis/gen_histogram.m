function [histval,flops_CIJ] = gen_histogram( lattice_size, max_keep_states, do_plot_in )
% [histval,flops_CIJ] = gen_histogram( lattice_size, max_keep_states, [do_plot_in] )
%
% Generate histogram of work in CIJ cell
%
do_plot = 0;
if (nargin >= 3),
  do_plot = do_plot_in;
end;


target_up = round( lattice_size/2);
target_down = round( lattice_size/2);

left_size = round( lattice_size/2);
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
    fake_interaction_matrix = ones( npatches, npatches );
    [flops_total,flops_CIJ] = cal_CIJ_flops( npatches, ...
                                             left_patch_size, ...
                                             right_patch_size, ...
                                             fake_interaction_matrix );
    flops_total_with_fake_interaction_matrix = flops_total;

    [flops_total,flops_CIJ] = cal_CIJ_flops( npatches, ...
                                             left_patch_size, ...
                                             right_patch_size, ...
                                             interaction_matrix );

   flops_total_with_interaction_matrix = flops_total;
   disp(sprintf('flops_total with interaction_matrix %g', ...
                 flops_total_with_interaction_matrix));
   disp(sprintf('flops_total with fake interaction_matrix %g', ...
                 flops_total_with_fake_interaction_matrix));
   



   flops_CIJ = sort( flops_CIJ(:));
   csum_flops_CIJ = cumsum( flops_CIJ );
   total_flops = csum_flops_CIJ( length(csum_flops_CIJ) );

   % -----------------------------------------------------
   % note use max(1,flops_CIJ) to avoid problem with log10
   % -----------------------------------------------------
   log10_flops_CIJ = log10( max(1,flops_CIJ(:)) );
   imax = ceil(max(log10_flops_CIJ));
   imin = floor(min(log10_flops_CIJ));

   [histval,binval] = hist( log10_flops_CIJ, imin:imax );
   if (do_plot),
      subplot(1,2,1);
      hist( log10_flops_CIJ, imin:imax ); 
      xlabel('work (log10)');
      ylabel('number of cells');
      title('histogram of work in CIJ cells');

      subplot(1,2,2);
      plot( 1:length(csum_flops_CIJ),csum_flops_CIJ );
      xlabel('number of cells');
      ylabel('total work');
      title('cumulative distribution of work ');
  
   end;

   do_table = 1;
   if (do_table),
     n = length( flops_CIJ );

     p = 10;
     percent = 100*(total_flops - csum_flops_CIJ( round((1-p/100)*n) ))/ total_flops;
     disp(sprintf('top %d %% of cells account for %g %% of total work', p, percent ));

     p = 5;
     percent = 100*(total_flops - csum_flops_CIJ( round((1-p/100)*n) ))/ total_flops;
     disp(sprintf('top %d %% of cells account for %g %% of total work', p, percent ));
        
     p = 1;
     percent = 100*(total_flops - csum_flops_CIJ( round((1-p/100)*n) ))/ total_flops;
     disp(sprintf('top %d %% of cells account for %g %% of total work', p, percent ));

     p = 0.5;
     percent = 100*(total_flops - csum_flops_CIJ( round((1-p/100)*n) ))/ total_flops;
     disp(sprintf('top %g %% of cells account for %g %% of total work', p, percent ));

         
    end;
