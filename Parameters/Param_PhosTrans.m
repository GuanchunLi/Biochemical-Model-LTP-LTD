function [paras] = Param_PhosTrans(paras)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
load('phos_trans_mat.mat');
paras.pprate_nop_rev = pprate_nop_rev;
paras.pprate_nop = pprate_nop;
paras.pprate_p_rev = pprate_p_rev;
paras.pprate_p = pprate_p;
paras.sprate_rev = sprate_rev;
paras.sprate = sprate;
end

