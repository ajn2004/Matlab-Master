function cuda_compile(fname)
clc
ipaths = {['-I' fullfile('C:\ProgramData\NVIDIA Corporation\CUDA Samples\v8.0\common\inc')]};
mex(ipaths{:},fname)