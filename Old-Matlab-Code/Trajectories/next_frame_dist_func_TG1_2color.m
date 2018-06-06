function [q,data_all,xf_all,yf_all,framenum_all,a0_all,total_molecules,nrat,green_sum,red_sum]=next_frame_dist_func_TG1_2color(filename,microns_per_pixel,dmax,dfactor,dfactor2,mol_type)

%mol_type=1;

load(filename,'q','xf_all','yf_all','framenum_all','a0_all','total_molecules','nrat','green_sum','red_sum');
xf_all0=xf_all;
yf_all0=yf_all;
framenum_all0=framenum_all;
total_molecules0=total_molecules;
a0_all0=a0_all;
nrat0=nrat;
green_sum0=green_sum;
red_sum0=red_sum;

% select either the red or the green molecules
D_rat_min=0.00; %if using red/(red+green)
D_rat_max=0.99;
C_rat_min=0.99;
C_rat_max=1.0;

Ng_min=50;
Nr_min=50;

if mol_type==1
    % select green ones
    constraint=(nrat<D_rat_max & nrat>D_rat_min & green_sum>Ng_min & red_sum>Nr_min);
    selected=find(constraint==1);
end
if mol_type==2
    % select red ones
    constraint=(nrat<C_rat_max & nrat>C_rat_min & green_sum>Ng_min & red_sum>Nr_min);
    selected=find(constraint==1);
end

size(xf_all)
size(xf_all0)
max(max(selected))

xf_all=xf_all0(selected);
yf_all=yf_all0(selected);
framenum_all=framenum_all0(selected);
a0_all=a0_all0(selected);
nrat=nrat0(selected);
green_sum=green_sum0(selected);
red_sum=red_sum0(selected);
total_molecules=length(xf_all);


data_all=zeros(total_molecules,13);

xf_um=xf_all*microns_per_pixel;
yf_um=yf_all*microns_per_pixel;

xmean=mean(mean(xf_um));
ymean=mean(mean(yf_um));
%reassign coordinates for all molecules detected wrt center.
xf_um=xf_um-xmean;
yf_um=yf_um-ymean;

xmin=min(min(xf_um))-0.5;     % x and y limits of the graph to be drawn
xmax=max(max(xf_um))+0.5;
ymin=min(min(yf_um))-0.5;
ymax=max(max(yf_um))+0.5;

framenum_all0=framenum_all;     % save original frame numbers
first_frame_number=framenum_all0(1);  
framenum_all=framenum_all-first_frame_number+1; % force the first frame number to be 1

maxframe=max(max(framenum_all));
total_pairs=0;  
% figure

for frame=1:maxframe-1
    clear xc yc xn yn index_npc index_npn;   
    npc=0;
    npn=0;
    
    % find molecules in current ("c") frame and next ("n") frame
    i_c=find(framenum_all==frame);   %indices for molecules in current frame
    i_n=find(framenum_all==frame+1); %indices for molecules in next frame
                                           
    if min(size(i_c))>0 %make sure there are some in current
        for i=1:max(size(i_c))
           npc=npc+1;             % number of molecules in the current frame 
           xc(npc)=xf_um(i_c(i)); % coords of molecules in the current frame
           yc(npc)=yf_um(i_c(i));
           index_npc(npc)=i_c(i);
        end
    end
    if min(size(i_n))>0 %make sure there are some in next
        for i=1:max(size(i_n)) 
            npn=npn+1;             %number of molecules in the next frame       
            xn(npn)=xf_um(i_n(i)); %coords of molecules in the next frame
            yn(npn)=yf_um(i_n(i));
            index_npn(npn)=i_n(i);
        end
    end
    if npc>0 && npn>0 %if found some in current frame and some in next frame
%         clf
%         hold on
%         plot(xc,yc,'.b')
%         plot(xn,yn,'.r')
      
        %loop through molecules in the current frame
        for i=1:npc      
            %iniial seeding of variables 
            d1=1e9;  %nearest neighbor in next frame    
            n1=-1;   %index for nearest neighbor in next frame
            d2=1e9;  %next-nearest neighbor in next frame 
            n2=-1;   %index for next-nearest neighbor in next frame
            dsf=1e9; %closest distance in same frame
            nsf=-1;  %index for nearest neighbor in same frame

            %find the nearest neighbor in the same frame
            for j=1:npc    %scan through all particles in the given frame
                if i ~= j  %don't compare the molecule with itself
                    dx=xc(i)-xc(j);
                    dy=yc(i)-yc(j);
                    d=sqrt(dx*dx+dy*dy); %distance to the j-th particle
                    if d<dsf   %if closer than the closest one found so far, make this one the closest one
                        dsf=d; %make the new closest distance d, 
                        nsf=j; %and keep track of which one it was
                    end
                end
            end
            %now scan through all molecules in the next frame
            for j=1:npn          
                dx=xc(i)-xn(j); %find distance between ith molecule in current 
                dy=yc(i)-yn(j); %frame and all others in the next frame
                d=sqrt(dx*dx+dy*dy);
                if d<d1    %if necessary re-asign nearest and next-nearest
                    d2=d1;
                    n2=n1;      
                    d1=d;
                    n1=j;
                elseif d<d2
                    n2=j;
                    d2=d;
                end
            end
            if n1>0 && d1<dmax && d2>dfactor*dmax && dsf>dfactor2*dmax
                dx1=xn(n1)-xc(i); %calculate x and y displacements between the one in the current frame
                dy1=yn(n1)-yc(i); %and nearest in next frame
                
                % dmax is the maximum distance between a molecule in
                % the given frame and a molecule in the next frame to
                % be considered possibly the same molecule

                % also, the next nearest molecule can't be too close
                % (must be less than dfactor*dmax) so we don't confuse
                % it with the actual molecule in the next frame

                % also, the distance to the nearest one in the same frame 
                % has to be large enough (dfactor2*dmax) so IT isn't confused 
                % with the molecule in the current frame

                % in principle, these constraints should minimize errors

                % if it passes all these constraints, save the pair
                % molecule is likely the same

%                 plot([xc(i) xn(n1)],[yc(i) yn(n1)],'g','LineWidth',4)   

                total_pairs=total_pairs+1;
                %data_all=[index_c index_n frame_c frame_n xc yc NN_dx NN_dy NN_dist_c NN_dist_n NNN_dist_n a0_c a0_n]
                data_all(total_pairs,1)=index_npc(i);
                data_all(total_pairs,2)=index_npn(n1);
                data_all(total_pairs,3)=frame+first_frame_number-1;
                data_all(total_pairs,4)=frame+first_frame_number;
                data_all(total_pairs,5)=xc(i)+xmean;
                data_all(total_pairs,6)=yc(i)+ymean;
                data_all(total_pairs,7)=dx1;
                data_all(total_pairs,8)=dy1;
                data_all(total_pairs,9)=dsf;
                data_all(total_pairs,10)=d1;
                data_all(total_pairs,11)=d2;
                data_all(total_pairs,12)=a0_all(index_npc(i));
                data_all(total_pairs,13)=a0_all(index_npn(n1));
            end
        end
%         axis equal
%         xlim([xmin xmax]);
%         ylim([ymin ymax]);
%         title(['Frames: ' num2str(frame) ' -' num2str(frame+1)]);
%         drawnow
%         hold off
    end
end

data_all=data_all(1:total_pairs,:);