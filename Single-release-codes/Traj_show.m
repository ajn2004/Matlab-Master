% Traj Show
% A script to display trajectories

% build trajectories

% Tolerance the intial  localizations
% ind = iters < 45 & N < fluor & xf_crlb.^0.5*q < 1 & yf_crlb.^0.5*q < 1;
ind = xf_crlb.^0.5*q < 1 & yf_crlb.^0.5*q < 1;
ind = N > 25& ind & zf_all*q < 0.45 & zf_all*q > -0.7;
% Tolerance logic, less iterations mean fit converged faster, indicating
% better data for fit model. Number of photons detected should be less than
% the total number of photons on the frame. Exclude data with uncertainty
% which is above a micron in X and Y
f = figure;
tg = uitabgroup(f);
% Build trajecs
xs = xf_all(ind);
ys = yf_all(ind);
zs = zf_all(ind);
zs = zs - min(zs);
fs = framenum_all(ind);
ind2 = find(fs/2 ~= round(fs/2)); % grab all odd number frames
ind3 = find(fs/2 == round(fs/2)); % grab all even number frames
[m,n,o] = size(sdi1);
t5 = uitab(tg,'Title','Frames');
t55 = uitabgroup(t5);
for i = 1:o
    ax = axes(uitab(t55,'Title',['F ', num2str(i)]));
    imagesc(ax,sdi1(:,:,i));
    hold on
    wind = -3:3;
    colormap(ax,'gray')
    [xf,yf,sx,sy,Neat,O] = mle_Gauss(sdi1((m+1)/2+wind,(n+1)/2+wind,i),110*pi/180);
    xfm(i) = xf + cents(i,1);
    yfm(i) = yf + cents(i,2);
        plot(ax,xf_all(i) - cents(i,1),yf_all(i)- cents(i,2),'bx');
    plot(ax,xf+n/2,yf+m/2,'rx');
    title(['Sx = ', num2str(sx),' Sy = ', num2str(sy), 'N = ', num2str(Neat),' O = ', num2str(O)]);
    hold off
    axis image
end
t2 = uitab(tg,'Title','Uncertainties');
t22 = uitabgroup(t2);
tx = uitab(t22,'Title','X and Y uncertainty');
ax = axes(tx);
histogram(ax,xf_crlb(ind).^0.5*q*1000);
hold on
histogram(ax,yf_crlb(ind).^0.5*q*1000);
hold off
xlabel('Uncertainty in nm')
ylabel('Frequency');
legend('Xf unc','Yf unc');
tz = uitab(t22,'Title','Axial uncertainty');
ax = axes(tz);
histogram(ax,zf_crlb(ind).^0.5*q*1000);
xlabel('Uncertainty in nm')
ylabel('Frequency');
t3 = uitab(tg,'Title','Number of Photons');
ax = axes(t3);
histogram(ax,N(ind));
xlabel('Photons Fitted');
ylabel('Frequency')
data = rat_view(sdi1);
t4 = uitab(tg,'Title','Traces and Processing');
ax = axes(t4);
plot(ax,data(end-1,:)./data(end,:));
hold on
plot(ax,fs,data(end-1,fs)./data(end,fs),'o');
hold off
xlabel('Framenumber');
ylabel('Local Sum / remaining sum')


t1 = uitab(tg,'Title','3D projection with Trajs');
ax = axes(t1);
set_scale_ax(std(sdi1,1,3),q,4,ax);
hold on
plot3(ax,xs(ind2),ys(ind2),zs(ind2),'r.');
plot3(ax,xs(ind3),ys(ind3),zs(ind3),'g.');
for i = 1:numel(ind2)
    chk = find(fs == fs(ind2(i)) + 1);
    if ~isempty(chk)
        dx = [xs(ind2(i)), xs(chk)];
        dy = [ys(ind2(i)), ys(chk)];
        dz = [zs(ind2(i)), zs(chk)];
        plot3(ax,dx,dy,dz,'k');
    end
end
axis image
ztix = floor(max(zs)/4);
set(gca,'Ztick',0:ztix:max(zs));
set(gca,'ZtickLabel',(0:ztix:max(zs))*q);
title('Image of Trajectories')
xlabel('Position in um')
ylabel('Position in um')
zlabel('Axial Position um')
hold off