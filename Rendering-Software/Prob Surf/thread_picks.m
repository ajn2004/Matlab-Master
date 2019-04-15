% Thread Picker
tots = 0;
for i = 7:15
    for j = 7:15
        for k = 7:15
            if i*j*k < 1024 && i*j*k / 32 == round(i*j*k/32)
                disp(['[',num2str(i),',',num2str(j),',',num2str(k),',] yields ', num2str(i*j*k), ' total threads'])
                if tots < (i*j*k) 
                    tops = [i,j,k];
                    tots = i*j*k;
                end
            end
        end
    end
end