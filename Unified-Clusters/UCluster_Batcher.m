%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unified Cluster Batcher
%
% This is a program to repeatedly run func_you_clusters
% The user will choose whether to loop over parameters or files
%
% AJN 8/25/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc;

perc_nn = 0.8;
min_num = 10;
% loop = 'params';  %if you want to loop parameters uncomment this line
loop = 'files';   % if you want to loop files uncomment this line
% loop = 'both';    %if you want to loop both parameters and files, uncomment this line.... you crazy cowboy




[fname, fpath] = uigetfile('*NG*');
cd(fpath);

if numel(perc_nn) > 1 || numel(min_num) >1
    load(fname,'data','w');
    cluster_id = [];
    count = 0;
    for i = perc_nn
        cluster_id_temp = [];
        for j = min_num
            count = count +1;
            [cluster_id_tempy] = func_uclusters(w, data(:,1),data(:,2),i, j);
            cluster_id_temp = [cluster_id_temp,cluster_id_tempy];
            count/(numel(perc_nn)*numel(min_num))
        end
        cluster_id = cat(3, cluster_id,cluster_id_temp);
    end
    finvar = func_scree(data(:,1),data(:,2),cluster_id);
    figure
    for i = 1:numel(min_num)
        subplot(1,2,1); plot(perc_nn,finvar(:,i,1));xlabel('perc_nn in um'); ylabel('variance in x');
        hold on
        subplot(1,2,2); plot(perc_nn,finvar(:,i,2));xlabel('perc_nn in um'); ylabel('variance in y');
        hold on
    end
    
    hold off
elseif numel(perc_nn) == 1 && numel(min_num) == 1
    finfo = dir('*NG*');
    figure('Units','Normalized','OuterPosition',[0 0 1 1]);
    for i = 1:numel(finfo)
        %         load(finfo(i).name,'xf_all','yf_all','q','w');
        load(finfo(i).name);
        xf_all = data(:,1);
        yf_all = data(:,2);
        q = 1;
        while true
            [cluster_id] = func_uclusters(w, data(:,1),data(:,2), perc_nn, min_num);
            scatter(xf_all(cluster_id >0),yf_all(cluster_id >0),10,cluster_id(cluster_id >0));
            colormap(lines(1000));
            axis image
            getout = 0;
            while true  %aggresive line of questioning
                disp(['Perc_nn was ', num2str(perc_nn),' and Min num was ', num2str(min_num)]);
                answer = input('Are these settings ok?(y/n) :','s');
                if strcmp(answer,'N') || strcmp(answer,'n')
                    ans2 = input('What was the problem, (P)erc_nn or (M)in_num? (P/M)?:', 's');
                    if strcmp(ans2,'p') || strcmp(ans2,'P')
                        perc_nn = input('Enter a new percentage(enter 0.1 fopr 10%) to scale the mean NN by: ');
                        break
                    elseif strcmp(ans2,'m' ) || strcmp(ans2,'M')
                        min_num = input('Enter a new minimum number of overlapping molecules: ');
                        break
                    else
                        msgbox('THIS IS NOT A GAME ANSWER THE QUESTION CORRECTLY!');
                    end
                elseif strcmp(answer,'Y') || strcmp(answer,'y')
                    save([finfo(i).name(1:end-4),'_clust.mat'],'data','w','cluster_id','perc_nn','min_num')
                    getout = 1;
                    break
                else
                    msgbox('ANSWER THE QUESTION CORRECTLY OR I WILL SHUT DOWN THE COMPUTER!')
                end
            end
            if getout == 1
                break
            end
        end
    end
end