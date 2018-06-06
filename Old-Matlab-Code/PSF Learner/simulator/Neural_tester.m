%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Neural Tester
%
% Tests the neural net results on frames of simulated data
%
% AJN 10/23/15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc

%% Declare variables
frames = 20;
max_mol = 20;
max_pho = 1000;
min_pho = 100;
bkg = 2;
pixx = 50;
pixy = 50;
q = 132;
NA = 1.45;
lam = 573;

load('Neural_thetas_3.mat');
% load('Neural_theta.mat');
idx = [];
mols = [];
mol_locs = [];
Ns = [];
h = waitbar(0,'0% Frames Created');
%% Build image by calling func_frame_sim.m
for i = 1:frames
    [i4(:,:,i), index, num_mol, mol_pix,N] = func_frame_sim(max_mol, max_pho, bkg, pixx, pixy, q, NA, lam, min_pho);
    waitbar(i/frames, h,[num2str(100*i/frames),'% Frames created']);
    idx = [idx; index];
    mols = [mols; num_mol];
    mol_locs = [mol_locs; mol_pix];
    Ns = [Ns; N(:)];
    waitbar(i/frames, h,[num2str(100*i/frames),'% Frames created']); 
end

close(h)

%Elements of rolling ball subrtraction copied from Einzel
rball=6; %radius of rolling ball
se = strel('ball',rball,rball,0); %structural element, i.e. rolling ball
FWHM=1; %FWHM of gaussian smoothing in pixels
rk=(FWHM)/sqrt(2*log(2)); %1/e^2 smoothing radius in pixels
kw=20; %kernal width of smoothing function
[Xgs,Ygs]=meshgrid(-kw/2:kw/2,-kw/2:kw/2);
kd=sqrt(Xgs.*Xgs+Ygs.*Ygs);
gs=exp(-2*kd.*kd/(rk*rk));
gs=gs/sum(sum(gs)); %smoothing function normalized to have area = 1
rbox = 3;
box_overlap_factor = 1.5; %if center to center closer, don't include either
w_mask = round(rbox*box_overlap_factor); %width of masking when a "high" pixel is identified

%% pad i4
% i5 = [zeros(pixy,3), i4, zeros(pixy,3)];
% i5 = [zeros(3,pixx+6); i5; zeros(3,pixx+6)];

locs = [];
num_o_mol = 0;

for i = 1:frames
%% Background subtract
%     i4_gs = uint16(conv2(i4(:,:,i),gs,'same')); %smoothed original
%     bkg = double(imopen(i4_gs,se));
%     
%     iprod=double(i4(:,:,i))-bkg;
%     i5=iprod.*(iprod>0); %set negative values to 0
    
    [loc, num_mols, a3] = neural_calc(i4(:,:,i), Theta1, Theta2);
    locs = [locs; loc];
    num_o_mol(i,1) = num_mols;
    a(:,:,i) = a3;
end
% pos(:,1) = locs(:,1) - 3;
% pos(:,2) = locs(:,2) - 3;

mols = [0; mols];
num_o_mol = [ 0; num_o_mol];
f1 = figure;
for i = 1:frames
    imagesc(i4(:,:,i));
    colormap(gray);
    hold on
    % plot(pos(:,1),pos(:,2),'.r', 'MarkerSize',15);
    plot(locs(sum(num_o_mol(1:i))+1:sum(num_o_mol(1:i+1)),1),locs(sum(num_o_mol(1:i))+1:sum(num_o_mol(1:i+1)),2),'.r', 'MarkerSize',15);
    plot(mol_locs(sum(mols(1:i))+1:sum(mols(1:i+1)),1),mol_locs(sum(mols(1:i))+1:sum(mols(1:i+1)),2),'.b','MarkerSize',15,'Marker','x');
    hold off
    drawnow
%     waitforbuttonpress
    M(i) = getframe(gcf, [55,45,500,400]);
end
movie(M)
% figure
% imagesc(a3);
% colormap(gray);
% title('a3');
% hold on
% % plot(pos(:,1),pos(:,2),'.r', 'MarkerSize',15);
% plot(locs(:,1),locs(:,2),'.r', 'MarkerSize',15);
% plot(mol_locs(:,1),mol_locs(:,2),'.b','MarkerSize',15,'Marker','x');
% hold off