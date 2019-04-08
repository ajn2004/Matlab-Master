
% function func_build_shift_correct(zin,frames,y,cyc)

% Cycle over offset to minimize standard deviation of height corrected
% positions. THis is an analytical equivalent of a sum of least squares
% calculation.
frames = framenumber;
zin = q*ncoords(:,3);
y = zavg;
for j = 1:cyc
    fnumber = frames +j;
for i = 1:numel(zin)
    correc(i) = y(mod(fnumber(i),cyc)+1);
    zout(i) = zin(i) - correc(i);

end
    soff(j) = std(zout);
    plot(frames,zout,'.')
    ylim(2*[min(zin),max(zin)])
    title(['Correction with Offset set to ', num2str(j*360/cyc),' degrees']);
    drawnow
    M(j) = getframe(gcf);
end
movie2gif(M,'C:\Users\AJN Lab\Dropbox\Presentations\4-5-19 Wrap Up (3d Shift Storm)\offset.gif','DelayTime',0.04,'LoopCount',Inf)
clear M
offs = find(soff == min(soff)); % grab the minimum value
sizes = 0.7:0.01:1.3;
fnumber = frames + offs;
for j = 1:numel(sizes)
    y1 = y*sizes(j);
    for i = 1:numel(zin)
    correc(i) = y1(mod(fnumber(i),cyc)+1);
    zout(i) = zin(i) - correc(i);

    end
    plot(frames,zout,'.')
    ylim(2*[min(zin),max(zin)])
    title(['Scale set to ',num2str(sizes(j)),'x'])
    drawnow
    M(j) = getframe(gcf);
    sscal(j) = std(zout);
end
movie2gif(M,'C:\Users\AJN Lab\Dropbox\Presentations\4-5-19 Wrap Up (3d Shift Storm)\scale.gif','DelayTime',0.04,'LoopCount',Inf)
scale = find(sscal == min(sscal));
y = y*sizes(scale);

for i = 1:numel(zin)
    correc(i) = y1(mod(fnumber(i),cyc)+1);
    zout(i) = zin(i) - correc(i);
end

plot(frames,zout,'.')