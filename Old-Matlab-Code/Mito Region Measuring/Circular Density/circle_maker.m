function [im1] = circle_maker(x,y,radius, i1)
% this function will take in 3 vectors of x, y, and r of equal size and
% plot circles of radius r(i) with centers (x(i),y(i))
% c will determine the color of the circle used

% AJN 9/22/15
[m,n] = size(i1);
im1 = i1;
xc = x + (n+1)/2;
yc = y + (m+1)/2;

% make the circle mathematically
for ii = xc-int16(radius):xc+(int16(radius))
    for jj = yc-int16(radius):yc+(int16(radius))
        tempR = sqrt((double(ii) - double(xc)).^2 + (double(jj) - double(yc)).^2);
        if(round(tempR) == double(int16(radius)))
            im1(ii,jj)=1;
        end
    end
end
% convert the circle to an image
% im1= i1;
end