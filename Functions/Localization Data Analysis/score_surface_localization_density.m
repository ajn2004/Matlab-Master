function score = score_surface_localization_density(data)
[Idx, D] = knnsearch(data,data,'k',3);
D = D(:,3);
mu_d = mean(D);
st_d = std(D);
clear dist mol_num
%Perform a neighbor search around each molecule using statistical
%parameters of nearest neighbor search
for i = 1:numel(data(:,1))    
    for j = 1:3
        dist(:,j) = (data(:,j) - data(i,j)).^2;
    end
    dist = (sum(dist,2)).^0.5;
    mol_num(i) = sum(dist < mu_d + 4*st_d);
end
score = exp(-(D-mu_d).^2/(2*st_d^2));
score = score.*(1-exp(-mol_num.'/mean(mol_num)));