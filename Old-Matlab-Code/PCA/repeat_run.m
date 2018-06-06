pix = 500;
cuda_compile('filtcut.cu')
i1 = ones(pix);
runs = 1000;
for i = 1:runs
    disp([num2str(100*i/runs),'% complete']);
    disp(['we have a mean value of ' num2str(mean(t)), ' and a SE of ',num2str(std(t)/numel(t)^0.5)]);
    [y,t(i)] = filter_cut(i1,3);
    drawnow
    clc
end
plot(t)
mean(t)