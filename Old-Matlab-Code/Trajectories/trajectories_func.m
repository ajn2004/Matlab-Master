load(matfile); %for batch script

c_index=dist_data(:,1);
n_index=dist_data(:,2);
xc=dist_data(:,5);
yc=dist_data(:,6);
dx=dist_data(:,7);
dy=dist_data(:,8);
big = size(c_index(:,1));
if big(1,1) > 1
npairs=max(max(size(xc)));
next_pair_in_chain=zeros(npairs,1);
n_connects=zeros(npairs,1);

for i=1:npairs %loop through current frames
    for j=1:npairs %loop through next frames
        if c_index(i)==n_index(j)
             n_connects(i)=n_connects(i)+1;
             n_connects(j)=n_connects(j)+1;
             next_pair_in_chain(j)=i;
             break
        end
    end
end

c_used=zeros(npairs,1);
n_traj=0;

for i=1:npairs
    if c_used(i)==0
        c_used(i)=1;
        n_traj=n_traj+1;
        trajectories(n_traj,1)=c_index(i);
        trajectories(n_traj,2)=n_index(i);
        next_i=next_pair_in_chain(i);
        column=3;
        while next_i ~= 0
            trajectories(n_traj,column)=n_index(next_i);
            c_used(next_i)=1;
            column=column+1;
            next_i=next_pair_in_chain(next_i);
        end
        traj_length(n_traj)=column-1;
    end
end
end
% figure
% set(gcf,'Name',[this_file(1:index(2)-1), ': trajectories']); %for batch
% plot(xc,yc,'.r'),axis image,hold on;
% for i=1:npairs
%     if n_connects(i)>=1
%       plot([xc(i) xc(i)+dx(i)],[yc(i) yc(i)+dy(i)],'Color',[0 1 0],'LineWidth',0.1);
%     else
%       plot([xc(i) xc(i)+dx(i)],[yc(i) yc(i)+dy(i)],'Color',[0 0 1]','LineWidth',0.1);
%     end
% end
% title('Trajectories');
% hold off

save(save_name);
if exist('trajectories','var')
    clear trajectories traj_length dist_data;
end
clear xf_all yf_all nrat red_sum green_sum;