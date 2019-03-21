function ind = start_w_z(zf, eps, z_cal)
[~,pam] = getdz(1,1,z_cal);
ind = zeros(numel(zf),1);
for i = 1:numel(zf)
    id = find(zf(i)-pam(:,1) == min(zf(i)-pam(:,1)));
    if eps(i) <=1.2+pam(id,2)/pam(id,3) && eps(i) >= pam(id,2)/pam(id,3)-0.2
        ind(i) = 1;
    else
        disp(eps(i));
    end
end