% close all
% figure('Units','Normalized','Outerposition',[0 0 1 1]);
count = 1;
clear t t2
for i = 0:0.001:10
ind = wtf > i;
ind2 = off_all > i;
% subplot(1,2,1); plot(xf_all(ind), yf_all(ind), '.b'); title('abs');xlabel(num2str(sum(ind)/numel(xf_all)));
% subplot(1,2,2); plot(xf_all(ind2), yf_all(ind2),' .b'); title('non abs');xlabel(num2str(sum(ind2)/numel(xf_all)));
t(count) = var(xf_all(ind));
x(count) = i;
% t2(count) = var(yf_all(ind));
drawnow;
count = count +1;
% waitforbuttonpress;
end

% figure
plot(x,t)
% hold on
% % plot(t2)
% hold off