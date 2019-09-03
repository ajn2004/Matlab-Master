function ajn_scatter(plts)
[m,n] = size(plts);

for i = 1:(n-1)
    plot(plts(:,1),plts(:,i+1),'.')
    hold on
end
hold off