function show_vector_field(ix,iy,iz)
% This function will show a vector field



[m,n,o] = size(ix);

for i = 1:m*n*o-1
    tic
    ii  = mod(i,m)+1;
    jj = mod(floor(i/m),n)+1;
    kk = floor(i/(m*n))+1;
    
    
    quiver3(ii,jj,kk,ix(ii,jj,kk),iy(ii,jj,kk),iz(ii,jj,kk));
    hold on
%     drawnow
    t(i) = toc;
ajn_wait(t,i,m*n*o);
end
axis equal
hold off