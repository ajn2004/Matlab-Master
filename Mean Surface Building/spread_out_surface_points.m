    function spread_surface = spread_out_surface_points(dense_surface, kk, m, iterations)
sdata = dense_surface;
sdata1 = sdata;

fixed = sdata(:,1) > max(sdata(:,1)) - 0.1 | sdata(:,1) < min(sdata(:,1)) + 0.1;
free = find(fixed == 0);
for i = 1:iterations
    [Idx, D] = knnsearch(sdata,sdata(free,:),'K',kk);
    for k = 1:numel(free)
        
        sub_set = sdata(Idx(k,2:end),:);
        r = sub_set - sdata(free(k),:);
        rh = -r./abs(r);
        coeff = pca(sub_set);
        proj_mat = coeff(:,3).*coeff(:,3).';
        orth_proj = (proj_mat*rh.').';
        rallowed = rh - orth_proj;
        kick = m*sum(rallowed);
        if sum(kick.*kick).^0.5 < 0.5
            sdata1(free(k),:) = sdata(free(k),:) + kick;
%             disp('kick')
        end
            
        
    end
    sdata = sdata1;
    
end
spread_surface = sdata;