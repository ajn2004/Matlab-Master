p = 1;
clip = 2;
sumy =zeros(psize*2+1, psize*2+1);
for p = 1:numel(points(:,1))
ysub = i1(points(p,2) - psize : points(p,2) +psize,points(p,1) - psize : points(p,1) +psize,:);
sumy = sumy + ysub;
end

exp = 0.04;

x = 1/exp:1/exp:numel(ysub(1,1,:))/exp;
F = x/numel(ysub(1,1,:));
% si1 = sum(sum(ysub));
si1 = sum(sum(sumy))/p;
si1 = si1(:);
Fi1 = abs(fft(si1)/numel(x));
F1 = (Fi1(1:round(numel(si1)/2)));
F1(2:end) = 2*(F1(2:end));
plot(F(clip:end - round(numel(si1)/2)),F1(clip:end))
% plot(Pi1(1:round(end/2)))

