% ---------------------------------
% simple script to explore work on patches and
% potential load imbalance
%
% assume  equal split of 
% number of left sites and number of right sites
% ---------------------------------
nthreads = 16;
nsites = 32;
nsites_left = nsites/2;
nsites_right = nsites - nsites_left;

ntarget_up = 6;
ntarget_down = 6;
keptstates = 4000;

states_left = zeros( (ntarget_up+1), (ntarget_down+1));
states_right = zeros( (ntarget_up+1), (ntarget_down+1));


for nleft_down=0:ntarget_down,
for nleft_up=0:ntarget_up,
   nright_up = ntarget_up - nleft_up;
   nright_down = ntarget_down - nleft_down;

   states_left(nleft_up+1,nleft_down+1) = nchoosek( nsites_left, nleft_up)*nchoosek( nsites_left, nleft_down); 

   states_right(nright_up+1,nright_down+1) = nchoosek( nsites_right, nright_up)*nchoosek( nsites_right, nright_down);
end;
end;

states_left = states_left/sum( sum(states_left));
states_left = keptstates * states_left;
states_left = ceil( states_left );

states_right = states_right/sum( sum(states_right));
states_right = keptstates * states_right;
states_right = ceil( states_right );


figure(1);
clf
mesh( states_left ); 
title(sprintf('left states,ntarget-up=%d,nsites=%d',ntarget_up,nsites));
print -dpng fig1.png


patchsizes = reshape( states_left, [prod(size(states_left)),1]);
npatches = (ntarget_up+1)*(ntarget_down+1);

work = zeros(npatches,npatches);
work_diag = zeros(npatches,1);

ipatch = 0;
for nleft_down=0:ntarget_down,
for nleft_up=0:ntarget_up,

   nright_down = ntarget_down - nleft_down;
   nright_up = ntarget_up - nleft_up;

   ipatch = 1 + (nleft_up + nleft_down * (ntarget_up+1));

   for nleft_down2=0:ntarget_down,
   for nleft_up2=0:ntarget_up
     nright_down2 = ntarget_down - nleft_down2;
     nright_up2 = ntarget_up - nleft_up2;

     jpatch = 1 + (nleft_up2 + nleft_down2 * (ntarget_up+1));

     nrowA = states_left( nleft_up+1,nleft_down+1);
     ncolA = states_left( nleft_up2+1,nleft_down2+1);

     nrowB = states_right( nright_up+1,nright_down+1);
     ncolB = states_right( nright_up2+1,nright_down2+1);
     % -------------------
     % Ak = nrowA by ncolA
     % Bk = nrowB by ncolB
     % -------------------

     if (ipatch == jpatch),
       % --------------------------
       % kron( HL, eye)*X = X * HL'
       % kron( eye, HR)*X = HR * X
       % --------------------------
       flops1 = nrowA * ncolA * ncolB;
       flops2 = nrowB * ncolB * ncolA;
       flops = flops1  + flops2; 
       work(ipatch,jpatch) = flops;
     else

       flops1 = ncolA * nrowB * ( nrowA + ncolB );
       flops2 = nrowA * ncolB * ( ncolA + nrowB );
       flops = min(flops1,flops2);
       work(ipatch,ipatch) = work(ipatch,ipatch) + flops;
     end;

   end;
   end;
end;
end;

figure(2);
clf;
plot( 1:npatches, sort(sum(work)), 'bx-'); 
title(sprintf('estimated work, ntarget=%d,nsites=%d ',ntarget_up,nsites));
print -dpng fig2.png


% distribute work
weights = sort( sum(work) );
thread_load = zeros(nthreads,1);

   for ipatch=1:npatches,
      % find the thread with least load
      [min_val, min_thread] = min( thread_load);
      thread_load(min_thread) = thread_load(min_thread) + weights(npatches-ipatch+1);
   end;

figure(3);
clf
avg_load = sum(thread_load)/nthreads;
plot( 1:nthreads, sort(thread_load),'*', 1:nthreads, avg_load * ones(nthreads,1), 'r-');
title(sprintf('redistributed work load,nthreads=%d',nthreads));
print -dpng fig3.png

   
    
