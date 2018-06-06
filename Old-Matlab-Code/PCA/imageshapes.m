% Image reshaping
w = 3;
[m,n] = size(i1);

ms = floor(m/(2*w+1));
ns = floor(n/(2*w+1));
x = [];
tic
for i = w+1:(2*w+1):m
    for j = w+1:(2*w+1):n
        i2 = i1(i-w:i+w,j-w:j+w);
        x = [x;i2(:).'];
    end
end
toc