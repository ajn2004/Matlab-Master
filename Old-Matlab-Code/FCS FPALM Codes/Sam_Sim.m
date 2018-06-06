%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation for Sam
%
% This is a simulation meant to compare rotational dynamics of a particle
% free to rotate and one that is not
%
%
% AJN 1/4/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc

%%%%%%%%%%%%%%% From Sam Hess %%%%%%%%%%%%%%%%%%%%%%%%%%%
np=1000000; 
dtheta=0.001; 
thetaval=0:dtheta:pi; 
ptheta=sin(thetaval); 
intptheta=0.5*(1-cos(thetaval)); 
plot(thetaval,intptheta); 
  
for i=1:np 
    s=rand; 
    ind=find(intptheta>s,1,'first'); 
    theta(i)=thetaval(ind); 
end 
  
%theta=rand(np,1)*pi; 
  
figure
phi=rand(np,1)*2*pi; 
phi = phi.';
xp=sin(theta).*cos(phi); 
yp=sin(theta).*sin(phi); 
zp=cos(theta); 
plot3(xp,yp,zp,'.b','MarkerSize',0.1); 
axis equal

%%%%%% End Sam Part %%%%%%%%%%%%%%%%


full = mean(zp.*zp);

for i = 1:101
parts(i) = mean(zp(cos(theta).^2 > (i-1)/100).^2);
end
plot(0:0.01:1,parts/full)