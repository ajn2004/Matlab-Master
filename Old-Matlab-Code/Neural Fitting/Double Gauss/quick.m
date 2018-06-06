ind = find(abs(results) > 0.5);
for i = 1:numel(ind)
    distances(i) = ((truth(ind(i),2) - truth(ind(i),4))^2 + (truth(ind(i),3) - truth(ind(i),5))^2);

end
histogram(distances*133)