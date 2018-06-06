% Temporal Test for PCA

pad = 3;
% pix = 500;
count = 1;
x = 0:0.01:7;
ks = 10.^x;
for k = ks
    100* count / numel(ks)
    i1 = ones(round(k));
    [m,n] = size(i1);
    i2 = [zeros(pad,n+2*pad);zeros(m,pad),i1,zeros(m,pad);zeros(pad,n+2*pad)];
    tic
    % image filter
    count = 1;
    for i = pad +1:m+pad
        for j = pad+1:n+pad
            i3 = i2(i-pad:i+pad,j-pad:j+pad);
            ims(count,:) = i3(:).';
            count = count +1;
        end
    end
    % end
    xl(count) = round(k);
    cpu(count) = toc; %end cpu version
    tic
    y = filter_cut(i1,pad);
    gpu(count) = toc;
    count = count +1;
end
hold on
plot(xl.^2*4,cpu,'.r');
plot(xl.^2*4,gpu,'.b');
legend('cpu','gpu');
hold off