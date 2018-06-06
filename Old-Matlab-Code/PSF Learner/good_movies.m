%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Good MOvies
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc
count = 1;

max_i = 100;
max_j = 100;
load('learned_theta.mat');
[xpix, ypix]= meshgrid(-3:3,-3:3);
x0true = 0;
y0true = 0;
N0 = 1000;
sigma2 = 2;   %this would be rayleigh radius in pixel space
sigma = sigma2/2;   % gaussian sigma
B = 3;
% w2 = waitbar(0, ' Creating points');
i1 = xpix.*0;
% Create a gaussian
figure 
hold on
cmap = ['y', 'm', 'c', 'r', 'g', 'b', 'w', 'k'];
for l = 1:8
    N = N0*(l/8);
for i = 1:max_i
%     x0true = -5 + 10*i/max_i;
%     y0true = x0true
    B = N/(i);
%     sigma = (i/max_i)*sigma2/2;
%     N = N0*(i/max_i);
    i1 = N.*(2*pi*(sigma)^2)^-1.*exp(-((xpix-x0true).^2 +(ypix-y0true).^2)./(2*(sigma).^2))+B^2;
    for j = 1:max_j
        i2 = imnoise(uint8(i1),'poisson');
        i3 = [ 1, double(i2(:)).'];
        prob(i,j)= sigmoid(i3*theta);
        
%         if prob(i,j) > 0.7
%                 imagesc(i2)
%                 colormap('gray')
%                 title(['Passed at ', num2str(prob(i,j)),'%']);
%                 xlabel(['Average value is ', num2str(i),' photons'])
%                 ylabel([num2str(100*j/max_j),'% of this value ', num2str(100*i*j/(max_i*max_j)), '% total']);
%                 drawnow
%         end
        disp([num2str(100*j/max_j),'% of this value ', num2str(100*((i-1)*max_j +j)/(max_i*max_j)), '% total'])
        count = count+1;
    end
    [h(i,:) , x] = hist(prob(i,:),0:0.01:1);
end

plot(1:max_i,mean(prob.'),cmap(l))
xlabel('Number of Photons');
ylabel('Sigmoid score');
end
figure
imagesc(i1);
colormap('gray');
axis image