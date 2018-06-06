cd('K:\Data\2-23-16 live mito tcs\Analysis\Toleranced\Centers\RDFS')
finfo = dir('*rdf.mat');
for i = 1:numel(finfo)
    load(finfo(i).name);
    subplot(5,5,i); plot(R,mean(g,2));
    title(finfo(i).name(1:end-4));
    xlabel('R in nanometers')
    ylabel('g(r)');
    xlim([0,1000])
end