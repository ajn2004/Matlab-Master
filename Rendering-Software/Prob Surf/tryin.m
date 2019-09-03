lines = 10;
fpath = 'C:\Users\AJN Lab\Dropbox\Lab Meetings AJN\May 1 2019';
close all
mols = 500;
x = rand(mols,1)*(lines-1) +1;
y = rand(mols,1)*(lines-1) +1;

c = plot(x,y,'.r')
axis equal
hold on
saveas(c, [fpath,'scatteredpnts.png'],'png');

for i = 1:lines
 c=  plot([1,lines],[i, i],'b')
hold on
end

for i = 1:lines
   c=plot([i, i],[1,lines],'b')
hold on
end

saveas(c, [fpath,'gridpoints.png'],'png');

c = plot(5.5,5.5,'bx','MarkerSize',10);

saveas(c, [fpath,'xmarksthegridpoints.png'],'png');

circ = 10000;
theta = 0:2*pi/circ:2*pi;
r = 1.3;
xc = r*cos(theta) + 5.5;
yc = r*sin(theta) + 5.5;
c = plot(xc,yc,'b')

saveas(c, [fpath,'circxmarksthegridpoints.png'],'png');
c = plot([5.5, 6.8], [5.5, 5.5],'b')
saveas(c, [fpath,'radcircxmarksthegridpoints.png'],'png');
