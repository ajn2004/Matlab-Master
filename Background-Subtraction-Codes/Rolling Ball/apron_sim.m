%% Sim for sep filter
% The goal is to have a script that will treat an image as a one
% dimensional array yet adequately build an apron along the edges by using
% the modulo function to detect when it is nearby
m=100;
n = 123;
o = 1;
im1 = ones(m,n,o);
 pixw = 5;
 
 for i = 0:numel(im1(:))+pixw
   
 end
 im2 = reshape(i2,m+2*pixw+1,n);