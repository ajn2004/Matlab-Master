% K-Factor 1 D toy example
clearvars
clc
close all

x = 1:0.1:100;
y = 2*sin(x) + randn(1,numel(x))+5*exp(-(x-20).^2);
y = y + min(y)+50;
ny = y/max(y);
plot(x,ny)
M = 80;
%% Finding the minimum k is done by finding the lowest value of the image
my = min(ny);

%create a k space to search over
ks = 0.001:0.001:1-0.001;
ks = ks.';
k = 0.8;
% write up a p function
for j = 1:numel(ks)
    
    fy = ny;
    for i = 1:M
        
        % close all
        
        
        g(i,:) = fy >= (1+k^i)^-1;
        %     g(i,find(fy == my)) = 0;
        f(i,:) = (1 + k^i.*g(i,:))./(1+k^i);
        fy = fy./f(i,:);
    end
    res(j) = sum((ny - prod(f,1)).^2);
end

% figure
for i = 1:10
    subplot(10,2,1+(i-1)*2);plot(prod(f(1:2*i,:),1)); title('in');ylabel(num2str(2*i));
    subplot(10,2,2+(i-1)*2);plot(prod(f(2*i+1:end,:),1));title('out');
end


% determine k from minimum I
