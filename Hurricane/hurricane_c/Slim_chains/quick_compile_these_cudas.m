files = dir('*.cu')
for i = 2:numel(files)
    mexcuda(files(i).name)
end