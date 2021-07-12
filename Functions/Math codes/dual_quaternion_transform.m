function [R, T] = dual_quaternion_transform(data_0,data_1,w)
% This function will determine the rotational and transformational
% transform between data_0 and data_1 given the weighting vector w
if ~exist('w','var')
    w = data_0(1,:)*0+1;
end
C1 = zeros(4);
C2 = zeros(4);
Ident = [1 0 0; 0 1 0; 0 0 1];
% Build final matrix
for i = 1:numel(data_0(1,:))
    xi = data_1(:,i);
    yi = data_0(:,i);
    kx = [ 0 -xi(3) xi(2); xi(3) 0 -xi(1); -xi(2) xi(1) 0];
    ky = [ 0 -yi(3) yi(2); yi(3) 0 -yi(1); -yi(2) yi(1) 0];
    C1 = C1 + w(i)*[ky*kx+yi*xi.', -ky*xi; -yi.'*kx, yi.'*xi];
    C2 = C2 + w(i)*[-kx-ky xi-yi; -(xi-yi).' 0];
end
C2 = 2*C2;
C1 = -2*C1;
W = sum(w);
A = 0.5*((1/(2*W))*C2.'*C2 - C1 -C1.');
[V, D] = eig(A);
q = V(:,end);
qu = q(1:3);
kq = [ 0 -qu(3) qu(2); qu(3) 0 -qu(1); -qu(2) qu(1) 0];
R = (q(4)^2 - qu.'*qu)*Ident + 2*qu*qu.' + 2*q(4)*kq;
s = -C2*q/(2*W);
Wq = [q(4)*Ident - kq qu; -qu.' q(4)];
p = Wq.'*s;
T = p(1:3);
end

