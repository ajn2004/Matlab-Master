

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make a movie of an object that rotates in 3D
%
% AJN 10/28/15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oh = findobj(gca, 'type','surface');
clear M
count = 1;
for i = 0:3:360
    view(i*1,0);
    M(count) = getframe(gca, [-75 -60 600 580]);
    count = count +1;
end
for j = 0:3:360
    if j <91 || j > 269
    view(0,j);
    M(count) = getframe(gca, [-75 -60 600 580]);
    count = count +1;
    else
        view(180,j)
        M(count) = getframe(gca, [-75 -60 600 580]);
        count = count +1;
    end
end

for k = 0:3:360
    if k <91 || k > 269
        view(k,k)
        M(count) = getframe(gca, [-75 -60 600 580]);
        count = count +1;
    else
        view(k+180,k)
        M(count) = getframe(gca, [-75 -60 600 580]);
        count = count +1;
    end

end
movie2avi(M, 'rotated surface5.avi')