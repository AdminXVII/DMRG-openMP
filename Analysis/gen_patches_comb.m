function [left_patch_size,right_patch_size,interaction_matrix] =  ...
          gen_patches_comb( left_size, right_size, ...
		       target_up, target_down, ...
                       keep_left_states, keep_right_states)
% ------------------------------------------------
%  [left_patch_size,right_patch_size,interaction_matrix] =  ...
%          gen_patches_comb( left_size, right_size, ...
%		       target_up, target_down, ...
%                      keep_left_states, keep_right_states)
%
%  estimate the size of patches based on 
%  combinatorial arguments
% ------------------------------------------------

max_left_up    = min( left_size,  target_up);
max_left_down  = min( left_size,  target_down);
max_right_up   = min( right_size, target_up);
max_right_down = min( right_size, target_down);


% ---------------------------------------------------------
% setup logical mask to identify valid (nleft_up,nleft_down)
% and valid matching (nright_up,nright_down)
% ---------------------------------------------------------
isvalid_left = zeros(1+max_left_up,1+max_left_down);
isvalid_right = zeros(1+max_right_up,1+max_right_down);

for nleft_up=0:max_left_up,
for nleft_down=0:max_left_down,
 
    nright_up = target_up - nleft_up;
    nright_down = target_down - nleft_down;
 
    isvalid = (0 <= nright_up) && (nright_up <= max_right_up) && ...
              (0 <= nright_down) && (nright_down <= max_right_down);
    isvalid_left( 1 + nleft_up,1+nleft_down) = isvalid;
    if (isvalid),
      isvalid_right(1+nright_up,1+nright_down) = isvalid;
    end;
end;
end;

npatches = sum( isvalid_left(:) );

left_patch_size  = zeros(npatches,1);
right_patch_size = zeros(npatches,1);




left_ways  = zeros( max_left_up+1, max_left_down+1);
right_ways = zeros( max_right_up+1, max_right_down+1);

% ---------------------------------------------
% vectorized way to compute
% Cmat_left(1+nleft_up,1+nleft_down) = ...
%               nchoosek(left_size,nleft_up) *  ...
%               nchoosek(left_size,nleft_down)
% ---------------------------------------------
nleft_up     = 0:max_left_up;
nleft_down   = 0:max_left_down;
C_nleft_up   = bincoeff( left_size, nleft_up);
C_nleft_down = bincoeff( left_size, nleft_down);

%   ------------------------------------------------------
%   Cmat_left = C_nleft_up(:) * transpose(C_nleft_down(:));
%   ------------------------------------------------------

for nleft_down=0:max_left_down,
for nleft_up=0:max_left_up,
   if (isvalid_left(1+nleft_up,1+nleft_down)),
%    ---------------------------------
%    note, beware of possible overflow
%    or loss of accuracy
%    since nchoosek(n,k) can quickly become
%    very large numbers
%
%    for example, nchoosek(144,72) = 1.4802*10^42
%    ---------------------------------
     left_ways(1+nleft_up,1+nleft_down) =  ...
         C_nleft_up(1+nleft_up) * ...
         C_nleft_down(1+nleft_down);

   end;
end;
end;



% ---------------------------------------------
% vectorized way to compute
% Cmat_right(1+nright_up,1+nright_down) = ...
%               nchoosek(right_size,nright_up) *  ...
%               nchoosek(right_size,nright_down)
% ---------------------------------------------
nright_up     = 0:max_right_up;
nright_down   = 0:max_right_down;
C_nright_up   = bincoeff( right_size, nright_up);
C_nright_down = bincoeff( right_size, nright_down);

%   ----------------------------------------------------------
%   Cmat_right = C_nright_up(:) * transpose(C_nright_down(:));
%   ----------------------------------------------------------

for nright_down=0:max_right_down,
for nright_up = 0:max_right_up;
   if (isvalid_right(1+nright_up,1+nright_down)),
     right_ways( 1+nright_up,1+nright_down) = ...
         C_nright_up(1+nright_up) * ...
         C_nright_down(1+nright_down);
   end;
end;
end;
     




total_left_ways = sum( left_ways(:) );
total_right_ways = sum( right_ways(:) );


% --------------------------------------------------------------------------------------
%
% set use_ceil = 1  to generate more patches
% set use_ceil = 0  to generate fewer patches
%
% --------------------------------------------------------------------------------------
% use ceiling() function  would generate at least 1 row per valid (nleft_up,nleft_down)
% this may generate more patches
%
% use round() function  would remove many unlikely or small states
% this may generate fewer patches
%
% using more patches would spread the rows across more patches and
% may perform less work overall
% --------------------------------------------------------------------------------------
use_ceil = 0;
if (use_ceil),
  left_ways = ceil((left_ways/total_left_ways) * keep_left_states );
  right_ways = ceil((right_ways/total_right_ways) * keep_right_states);
else
  left_ways = round((left_ways/total_left_ways) * keep_left_states );
  right_ways = round((right_ways/total_right_ways) * keep_right_states);
end;

nleft2patch = zeros( max_left_up+1,max_left_down+1);
nright2patch = zeros( max_right_up+1,max_right_down+1);

max_npatches = max( (max_left_down+1)*(max_left_up+1), ...
                    (max_right_down+1)*(max_right_down+1) );
patch2nleft_up = -1*ones( max_npatches,1);
patch2nleft_down = -1*ones( max_npatches,1);

npatches = 0;
for nleft_down=0:max_left_down,
for nleft_up=0:max_left_up,
  nright_up = target_up - nleft_up;
  nright_down = target_down - nleft_down;
  
  
  if (isvalid_left(1+nleft_up,1+nleft_down)),
   lsize = left_ways(1+nleft_up,1+nleft_down);
   rsize = right_ways(1+nright_up,1+nright_down);

   has_work = (lsize*rsize >= 1);
   if (has_work),
     npatches = npatches + 1;
     left_patch_size(npatches) = lsize;
     right_patch_size(npatches) = rsize;

     nleft2patch(1+nleft_up,1+nleft_down) = npatches;
     nright2patch(1+nright_up,1+nright_down) = npatches;

     patch2nleft_up(npatches) = nleft_up;
     patch2nleft_down(npatches) = nleft_down;
   end;

  end;

end;
end;



% -----------------------------------------------------------------
% Note, due to rounding and truncation, 
% sum(left_patch_size(:)) may not exactly equal to keep_left_states
% -----------------------------------------------------------------

need_adjust = 0;
if (need_adjust),
  

  % ----------------
  % adjust left part
  % ----------------
  while (sum(left_patch_size) ~=  keep_left_states),
     idiff = keep_left_states - sum(left_patch_size);

     if (idiff < 0),
        % -----------------------------
        % substract from largest states
        % -----------------------------
        [dummy,idx] = sort( left_patch_size,'descend');
     else
        % ----------------------
        % add to smallest states
        % ----------------------
        [dummy,idx] = sort( left_patch_size,'ascend');
     end;

     for i=1:npatches,
        if (idiff == 0),
            break;
        end;
  
        ipatch = idx(i);
        if (idiff  < 0),
             % ---------------------------
             % note avoid removing a patch
             % ---------------------------
             if (left_patch_size(ipatch) > 1),
               left_patch_size(ipatch) = left_patch_size(ipatch)-1;
               idiff = idiff + 1;
             end;
        else
           % -------------------------------
           % note avoid creating a new patch
           % -------------------------------
           if (left_patch_size(ipatch) >= 1),
             left_patch_size(ipatch) = left_patch_size(ipatch) + 1;
             idiff = idiff - 1;
           end;
        end; 
     end; % for

  end; % while

  isok = (sum(left_patch_size) == keep_left_states);
  if (~isok),
    error(sprintf('gen_patches_comb:sum(left_patch_size)=%d,keep_left_states=%d', ...
                               sum(left_patch_size),   keep_left_states));
  end;



  % ----------------
  % adjust right part
  % ----------------
  while (sum(right_patch_size) ~=  keep_right_states),
     idiff = keep_right_states - sum(right_patch_size);

     if (idiff < 0),
        % -----------------------------
        % substract from largest states
        % -----------------------------
        [dummy,idx] = sort( right_patch_size,'descend');
     else
        % ----------------------
        % add to smallest states
        % ----------------------
        [dummy,idx] = sort( right_patch_size,'ascend');
     end;

     for i=1:npatches,
        if (idiff == 0),
            break;
        end;
  
        ipatch = idx(i);
        if (idiff  < 0),
             % ---------------------------
             % note avoid removing a patch
             % ---------------------------
             if (right_patch_size(ipatch) > 1),
               right_patch_size(ipatch) = right_patch_size(ipatch)-1;
               idiff = idiff + 1;
             end;
        else
           % -------------------------------
           % note avoid creating a new patch
           % -------------------------------
           if (right_patch_size(ipatch) >= 1),
             right_patch_size(ipatch) = right_patch_size(ipatch) + 1;
             idiff = idiff - 1;
           end;
        end; 
     end; % for

  end; % while

  isok = (sum(right_patch_size) == keep_right_states);
  if (~isok),
    error(sprintf('gen_patches_comb:sum(right_patch_size)=%d,keep_right_states=%d', ...
                               sum(right_patch_size),   keep_right_states));
  end;
  


end;
         

% -------------------------
% keep only non-zero patches
% -------------------------
idx_nonzero = find( (left_patch_size > 0) & (right_patch_size > 0) );
left_patch_size = left_patch_size(idx_nonzero);
right_patch_size = right_patch_size(idx_nonzero);


patch2nleft_up = patch2nleft_up(idx_nonzero);
patch2nleft_down = patch2nleft_down(idx_nonzero);
% ---------------------------------------------------------------
% determine the sparse interaction matrix between (ipatch,jpatch)
% ---------------------------------------------------------------
interaction_matrix = zeros(npatches,npatches);
for jpatch=1:npatches,
 nleft_up_j = patch2nleft_up(jpatch);
 nleft_down_j = patch2nleft_down(jpatch);
 for ipatch=1:npatches,
   nleft_up_i = patch2nleft_up(ipatch);
   nleft_down_i = patch2nleft_down(ipatch);

   njump = abs(nleft_up_i-nleft_up_j) + ...
           abs(nleft_down_i-nleft_down_j);
   interaction_matrix(ipatch,jpatch) = (njump <= 1);
 end;
end;

