% Thread Picker
tots = 0;
for i = 2:16
    for j = 2:i
        for k = 2:i
            if i*j*k <= 1024 && i*j*k / 32 == round(i*j*k/32) 
                disp(['[',num2str(i),',',num2str(j),',',num2str(k),'] yields ', num2str(i*j*k), ' total threads'])
                if tots < (i*j*k) 
                    tops = [i,j,k];
                    tots = i*j*k;
                end
            end
        end
    end
end