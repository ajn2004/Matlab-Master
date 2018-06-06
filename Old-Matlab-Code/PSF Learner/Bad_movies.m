%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bad MOvies
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
count = 1;

max_i = 100;
max_j = 100;

for i = 1:max_i
    for j = 1:max_j
        i2 = imnoise(uint8(i.*ones(7)),'poisson');
        i3 = [ 1, double(i2(:)).'];
        prob(i,j)= sigmoid(i3*theta);
        
        % if prob(count) > 0.7
%         imagesc(i2)
%         colormap('gray')
%         title(['Passed at ', num2str(prob(i,j)),'%']);
%         xlabel(['Average value is ', num2str(i),' photons'])
%         ylabel([num2str(100*j/max_j),'% of this value ', num2str(100*i*j/(max_i*max_j)), '% total']);
%         drawnow
        % end
        disp([num2str(100*j/max_j),'% of this value ', num2str(100*i*j/(max_i*max_j)), '% total'])
        count = count+1;
    end
    [h(i,:) , x] = hist(prob(i,:),0:0.01:1);
end