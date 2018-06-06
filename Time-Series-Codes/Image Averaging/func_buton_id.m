function [points] = func_buton_id(ifin, baseline)
isub = ifin(:,:,baseline-10:baseline+10);
load('buton_thetas.mat');
% reshape

xt = [];

for i = 1:numel(isub(1,:,1))
    xt = [xt;isub(:,i,:)];
end
x = zeros(numel(xt(:,1,1)),numel(xt(1,1,:)));
for i = 1:numel(xt(:,1,1))
    x(i,:) = xt(i,1,:);
end

%calculate
X1 = [ ones(numel(x(:,1)),1), x];
z2 = X1*theta1.';
a2 = [ones(numel(z2(:,1)),1), sigmoid(z2)];
z3 = a2*theta2.';
a3 = round(sigmoid(z3));

% index conversion [i,j] <=> [(i-1)*numel(isub(1,:,1) + j]

% rebuild a3 into an image

i3 = [];
for i = 1:numel(isub(1,:,1)) 
        i3 = [i3,a3((i-1)*numel(isub(:,1,1)) + 1:(i)*numel(isub(:,1,1)))];
end

% cluster and return
count = 1;
% i4 = imerode
[xgrid, ygrid] = meshgrid(-3:3,-3:3);
z = exp(-(xgrid.*xgrid + ygrid.*ygrid)/(2*0.45));
zn = z / sum(z(:));
i4 = conv2(i3,zn,'same');
se = strel('ball', 4, 1, 2);
iprod = (i3 - imopen(i4,se));
% imagesc(round(i3 - iprod))
i5 = round(i3 - iprod);
s = regionprops(logical(i5), 'Centroid');
points = cat(1,s.Centroid);
% save('testing.mat','i3');
% background subtraction

end