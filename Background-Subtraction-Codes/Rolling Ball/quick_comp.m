% Quick test between cpu and gpu

i1 = readtiff('C:\Users\AJN Lab\Dropbox\Data\3-6-19 JEK Cells JF-646\sub\Bottom\HILO-cell4-5mw-autofocus_1.tif');

tic
icpu = roball(i1,6,6);
ct = toc;

tic
igpu = gpu_rball(i1);
gt = toc;