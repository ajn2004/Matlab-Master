function cuda_compile(fname)
clc
ipaths = {['-I' fullfile('C:\ProgramData\NVIDIA Corporation\CUDA Samples\v8.0\common\inc')]};
if exist('mex_CUDA_win64.xml') ==2
    copyfile(['C:\Users\AJN Lab\Desktop\Code Development\Matlab Testing Folder\Neural Quhzx\Source Codes\Open this after you read the text file\mex_CUDA_win64.xml']);
end
mex(ipaths{:},fname)