function y = Update_stoc(y, idx, paras)
    i = idx(1); j = idx(2);
    k = paras.reaction_lst(i);
    y(k, j) = 1 - y(k, j);
end

