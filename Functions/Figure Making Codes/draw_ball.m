% Make a sphere
% We're going to make a 3 dimensional plot of points that form a theoretical 3Ball
% Inputs will be center, radius, number of points
% clearvars;
% center = [0 ,0, 1]; % Center in x,y
% radius = 2; % radius in normalized units
% point_number =  100000;
function draw_ball(center, radius, point_number, ax)
% Genereate random vectors, they'll point in random directions
for i = 1:point_number
    x(i) = rand-0.5;
    y(i) = rand-0.5;
    z(i) = rand-0.5;
    norm = (x(i)^2 + y(i)^2+ z(i)^2)^0.5;
    nx(i) = x(i)/norm;
    ny(i) = y(i)/norm;
    nz(i) = z(i)/norm;
end
sphere_x = nx*radius + center(1);
sphere_y = ny*radius + center(2);
sphere_z = nz*radius + center(3);
% plot3(nx,ny,nz,'.')
% axis equal
hold on
scatter3(ax, sphere_x,sphere_y,sphere_z,5,'filled')
hold off
end