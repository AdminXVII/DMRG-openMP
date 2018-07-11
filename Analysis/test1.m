% simple test 

do_plot = 0;
do_write = 1;
lattice_size = 64;
max_keep_states = 5000;

lattice_size_table(1) = 64; max_keep_states_table(1) = 5000;
%lattice_size_table(2) = 64; max_keep_states_table(2) = 2*5000;
%lattice_size_table(3) = 2*64; max_keep_states_table(3) = 5000;
%lattice_size_table(4) = 2*64; max_keep_states_table(4) = 2*5000;
%lattice_size_table(5) = 8*64; max_keep_states_table(5) = 5000;
%lattice_size_table(6) = 8*64; max_keep_states_table(6) = 2*5000;

ncase = length(lattice_size_table);

for icase=1:ncase,
   lattice_size = lattice_size_table(icase);
   max_keep_states = max_keep_states_table(icase);

   [flops_total_sweep] = sim_sweep( lattice_size, max_keep_states, do_plot, do_write);
   disp(sprintf('flops_total_sweep=%g,lattice_size=%d,max_keep_states=%d', ...
              flops_total_sweep,   lattice_size,   max_keep_states ));
end;
