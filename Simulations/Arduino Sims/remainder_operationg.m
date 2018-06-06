close all
clearvars
x = 4;
j = 1:4;
for i = 1:numel(j)
    x = x-1;
    x = rem(x,4);
    y(i) = x;
    x = x+4;
    if x < 0
        disp('IT HAPPENED');
        break;
    end
end
for i = 1:numel(j)
    x = x+1;
    x = rem(x,4);
    y(i+numel(j)) = x;
    x = x+4;
    if x < 0
        disp('IT HAPPENED');
        break;
    end
end
plot(1:2*numel(j),y)