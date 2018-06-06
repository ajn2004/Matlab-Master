clearvars -except times results hidden

for i = 1:numel(results(1,1,:))
    for j = 1:numel(results(1,:,1))
        thestds(i,j) = std(results(:,j,i));
    end
end
subplot(2,3,1); plot(hidden,thestds(:,1),'.b');title('Stdevs in xf- true');
subplot(2,3,2); plot(hidden,thestds(:,2),'.b');title('Stdevs in yf- true');
subplot(2,3,3); plot(hidden,thestds(:,3),'.b');title('Stdevs in N- true');
subplot(2,3,4); plot(hidden,thestds(:,4),'.b');title('Stdevs in sigx- true');
subplot(2,3,5); plot(hidden,thestds(:,5),'.b');title('Stdevs in sigy- true');
subplot(2,3,6); plot(hidden,thestds(:,6),'.b');title('Stdevs in b- true');

