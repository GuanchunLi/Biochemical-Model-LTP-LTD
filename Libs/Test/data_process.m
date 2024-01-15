load('phos_trans_mat_orig.mat');
A_normal = @(A) A - diag(sum(A, 1));
pprate_nop_rev = A_normal(triu(blue));
pprate_nop = A_normal(tril(blue));
pprate_p_rev = A_normal(triu(red));
pprate_p = A_normal(tril(red));
sprate_rev = A_normal(triu(black));
sprate = A_normal(tril(black));

save('phos_trans_mat.mat', 'pprate_nop', 'pprate_nop_rev', 'pprate_p',...
     'pprate_p_rev', 'sprate', 'sprate_rev');
