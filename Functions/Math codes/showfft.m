function [fx,ft] = showfft(x,t)
% a function to give the fast fourier transform and frequency space of
% variable array x and t respectively
fx = fft(x);
ft = fftime(t);
plot(fftshift(ft),abs(fftshift(fx)));
xlabel('Frequency (Hz)');
ylabel('AU');
title('FFT');
figure
plot(fftshift(ft),log(abs(fftshift(fx))));
title('Power Spectrum FFT')
ylabel('AU');
xlabel('Frequency (Hz)');