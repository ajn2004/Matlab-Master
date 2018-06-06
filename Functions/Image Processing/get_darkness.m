function get_darkness()
files = dir('*back*');
i1 = readtiff(files(1).name);
mi1 = mean(i1,3);
save('back_subtract.mat','mi1');
end