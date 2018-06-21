function draw_boxes(cents,pixw)
hold on

for i = 1: numel(cents(:,1))
    x = [cents(i,1) - pixw, cents(i,1) + pixw, cents(i,1) + pixw, cents(i,1) - pixw, cents(i,1) - pixw];
    y = [cents(i,2) - pixw, cents(i,2) - pixw, cents(i,2) + pixw, cents(i,2) + pixw, cents(i,2) - pixw];
    plot(x,y,'r');
end
hold off

