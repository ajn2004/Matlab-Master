function [finfo, fname, fpath]  = getfinfo(exts)
[fname, fpath] = uigetfile(['*.',exts]);
finfo = dir(['*.',exts]);
