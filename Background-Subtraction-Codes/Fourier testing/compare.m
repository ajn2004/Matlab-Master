count = 1;
for i = 100:1000
    x = rand(300,300,i);
    tic
    a = abs(cgpufourier1(x,1));
	t1(count) = toc;
    tic
    a = abs(cgpufourier2(x,1));
	t2(count) = toc;
    count = count +1;
end
plot(t1);
hold on
plot(t2);