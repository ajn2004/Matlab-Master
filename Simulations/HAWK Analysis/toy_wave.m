function i2 = toy_wave(i1)
[m,n,o] = size(i1);
H0 = [1 -1];
H1 = [1 1 -1 -1];
for i = 1:m
    for j = 1:n
        for k =1:2
            i2(i,j,k) = abs((i1(i,j,k)*H0(1) + i1(i,j,k+1)*H0(2)));
        end
        
    end
end