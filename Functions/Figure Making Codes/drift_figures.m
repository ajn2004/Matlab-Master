% Drift-Figures
close all
xf = (xf_all-mean(xf_all(1)))*q;
yf = (yf_all-mean(yf_all(1)))*q;
zf = (zf_all-mean(zf_all(1)))*q;
try
dx = (drifts(:,1) - mean(drifts(1,1)))*q;
dy = (drifts(:,2) - mean(drifts(1,2)))*q;
end
ff = framenum_all*0.02;
% subplot(1,2,1);
plot(ff,xf);
hold on
plot(ff,yf);
legend('X(t)','Y(t)')
xlabel('Time(s)');
ylabel('Position(um)');
title('X(t) and Y(t) Drift');
try
plot(ff,-dx)
plot(ff,-dy)
end
figure
[M,V] = eig(cov(xf,yf));
rots = M*[xf.';yf.'];
plot(ff,rots(2,:).')
title('Principle Component Position');
xlabel('Time (s)');
ylabel('Position (um)');
% plot(ff,zf-mean(zf));

% subplot(1,2,2);
figure
s = scatter(xf,yf,[],framenum_all*0.02);
axis equal
colormap('jet');
% s.MarkerFaceColor = s.MarkerEdgeColor;
xlabel('Position (um)');
ylabel('Position (um)');
title('X-Y Drift');
figure
plot(ff,zf)
xlabel('Time(s)');
ylabel('Position(um)');
title('Z(t) Drift');