% Example
% close all;
clearvars
clc
hold off
frames = 500;
pix = 7;
thresh = 0.1;

count =1;
N = 10;
rl =0.5
for N = 300:-10:0
    bkg = 5;
    
    [xgrid, ygrid] = (meshgrid(-(pix - 1)/2 : (pix-1)/2, -(pix - 1)/2: (pix - 1)/2));
%     for rl = 0:1
        % ims = zeros(frames,pix*pix);
        for i = 1:frames
            roll = rand;
            if roll > rl
                im = single(imnoise(uint16((N)*exp(-2*((xgrid - 2*(rand - 0.5) ).^2  + (ygrid - 2*(rand - 0.5) ).^2)./(3)) + bkg*rand),'Poisson'));
            else
                im = single(imnoise(uint16(bkg*rand(pix,pix)),'Poisson'));
            end
            ims(i,:) = im(:).';
        end
        
        for i = 1:pix*pix
            nims(:,i) = ims(:,i) - mean(ims(:,i));
        end
        
        % Get eigenvalues and vectors of covariance matrix
        [V, D] = eig(cov(nims));
        
        dmax = max(D(:));
        dthresh = dmax*thresh;
        [row, col] = find(D>dthresh);
        
        Vproj = V(:,col);
        
        vals = Vproj.'*nims.';
        
        % try
        %     scatter3(vals(1,:),vals(2,:),vals(3,:));
        % catch
            try
                histogram(vals(1,:).*vals(2,:));
            catch
                histogram(vals(1,:));
            end
        % end
        % hold on
        % end
        
        % for i = 1
        %     for j = i+1 : numel(vals(:,1))
        %
        %         plot(vals(i,:), vals(j,:),'.');
        hold on
        %     end
        % end
        % xlim([-25,25]);
        % ylim([-25,25]);
        title([' N is  ', num2str(N),'%']);
        % hold off
        drawnow
        M(count) = getframe(gcf);
        count = count +1;
    
%     end
end