function ft = fftime(t)
% This function will take a temporal array t(expected in seconds) and
% convert it into a frequency array ft(Hz)
dt = diff(t);
[m,n] = size(t);
ftt = (1/dt(1):1/dt(1):m/dt(1))/m;
% because the FFT is symmetric around the center we're going to try to
% match that here
sft = ftt(1:floor(m/2));
sft = sft(:);
if m/2 == round(m/2) % if m is even, the division should be easy
    ft = [sft;-sft(end:-1:1)]; % this should return a symmteric
else
    ft = [sft;ftt(floor(m/2)+1);-sft(end:-1:1)];
end
ft = reshape(ft,m,n);
