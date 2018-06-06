function T = clust_loop(im1, ranged,over)
% This function will create a Scree plot for clustering in im1


[row, col] = find(im1); % input im1 should be a logical image from an edge filter


Y = pdist([col,row]);
disp('got pdist');
Z = linkage(Y);
disp('Linkage');
count = 1;
for i = ranged
    tic
    T = cluster(Z,'cutoff',i,'criterion','distance');
    tib = toc;
    while tib < 0.1
        tib = toc;
    end
    scatter(col,-row,10,T,'Filled');
    colormap('Jet')
    title(['Cutoff is ', num2str(i)])
    axis image
    whitebg([0 0 0]);
    drawnow
    for j = 1:max(T)
        index = find(T==j);
        varx(j,count) = var(row(index));
        vary(j,count) = var(col(index));
    end
    finvar(count,1) = mean(varx(:,count));
    finvar(count,2) = mean(vary(:,count));
    count = count +1;
end

figure
plot(ranged,finvar(:,1),'.b')
hold on
plot(ranged,finvar(:,2),'.r')
hold off
if ~strcmp(over,'y')
while true
    ct = input('What Cutoff do you want? ');
    T = cluster(Z,'cutoff',ct,'criterion','distance');
    scatter(col,-row,10,T,'Filled');
    colormap('Jet')
    title(['Cutoff is ', num2str(i)])
    axis image
    whitebg([0 0 0]);
    drawnow
    s = input('Is this what you wanted? ','s');
        if strcmp(s,'Y') || strcmp(s,'y')
        break
        end
end
end

T = [T, row, col];