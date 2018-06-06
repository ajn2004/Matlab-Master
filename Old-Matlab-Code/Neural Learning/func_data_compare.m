function [missed_mol, false_p, victory, tot_num, N_fa, N_mi,N_p, off_p, off_mi, off_fa, sigs_mi  ] = func_data_compare(fname,truth_name, dist_thresh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Data examiner
%
% A script to analyze the symmetry between localized data and ground truth
% from simulations.
%
% AJN 5/2/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% User controleld variables
% truth_name = 'mock_data_pN_150_pS_1_bkg_4_mpf_5noverlap.mat';
% dist_thresh = 5; %this is the threshold distance in pixels between a data point and the simulated points to count as a succesful identification

%% User independent code
% [fname, fpath] = uigetfile('*.mat');

% cd(fpath);
load(truth_name);
load(fname);

% data(i,:) containes information for the ith molecule [coords(i,:), N(N_ind(i)), sigs(sigs_ind(i)), offset, frame, bkg];

frm_max = max(data(:,6));
victory = 0;
missed_mol = 0;
false_p = 0;
xfa = [];
yfa = [];
N_fa = [];
off_fa = [];
sigx_fa = [];
sigy_fa = [];
xmi = [];
ymi = [];
N_mi = [];
off_mi = [];
sigs_mi = [];
xpa = [];
ypa = [];
N_p = [];
off_p = [];
sigx_p = [];
sigy_p = [];
for i = 1:frm_max % loop over all frames of data
    indexr = find(framenum_all == i);
    indexs = find(data(:,6) == i);
    count = 0;
    % If the two indeces are non empty, loop over the simulated molecuels
    if numel(indexr) >0
        if numel(indexs) >0
            clear xf_sub yf_sub
            % create subset vectors
            xf_sub = xf_all(indexr);
            yf_sub = yf_all(indexr);
            N_sub = N(indexr);
            off_sub = off_all(indexr);
            sigx_sub = sigx_all(indexr);
            sigy_sub = sigy_all(indexr);
            
            % loop simulated molecules
            
            for j = 1: numel(indexs)
                dist_vec = ((xf_sub - data(indexs(j),1)).^2 + (yf_sub - data(indexs(j),2)).^2).^0.5;
                
                % counts the number of molecules within a distance from
                % true molecule
                for k = 1:numel(dist_vec)
                    if dist_vec(k) < dist_thresh
                        victory = victory + 1;
                        count = count +1;
                        xpa = [xpa; xf_sub(k)];
                        ypa = [ypa; yf_sub(k)];
                        N_p = [N_p; N_sub(k)];
                        off_p = [off_p; off_sub(k)];
                        sigx_p = [sigx_p ; sigx_sub(k)];
                        sigy_p = [sigy_p; sigy_sub(k)];
                        
                        sigy_sub(k) = -1;
                        xf_sub(k) = -1;
                        yf_sub(k) = -1;
                        N_sub(k) = -1;
                        off_sub(k) = -1;
                        sigx_sub(k) = -1;
                        sigy_sub(k) = -1;
                        indexs(j) = 0;
                    end
                end
            end
            false_p = false_p + (numel(xf_sub) - count); % count the false positives
            missed_mol = missed_mol + (numel(indexs)  - count);
            indexf = N_sub > -1;
            if sum(indexf) >0
                xfa = [xfa;xf_sub(indexf)];
                yfa = [yfa;yf_sub(indexf)];
                N_fa = [N_fa; N_sub(indexf)];
                off_fa = [off_fa ; off_sub(indexf)];
                sigx_fa = [sigx_fa; sigx_sub(indexf)];
                sigy_fa = [sigy_fa; sigy_sub(indexf)];
            end
            if sum(indexs >0) >0
                ind = find(indexs >0);
                for lk = 1:numel(ind)
                    xmi = [xmi;data((indexs(ind(lk))),1)];
                    ymi = [ymi;data((indexs(ind(lk))),2)];
                    N_mi = [N_mi; data((indexs(ind(lk))),3)];
                    off_mi = [off_mi ; data((indexs(ind(lk))),5)];
                    sigs_mi = [sigs_mi; data((indexs(ind(lk))),4)];
                    
                end
            end
        end
    else
        if sum(indexs) >0
            missed_mol = missed_mol + numel(indexs);
            ind = find(indexs >0);
                for lk = 1:numel(ind)
                    xmi = [xmi;data((indexs(ind(lk))),1)];
                    ymi = [ymi;data((indexs(ind(lk))),2)];
                    N_mi = [N_mi; data((indexs(ind(lk))),3)];
                    off_mi = [off_mi ; data((indexs(ind(lk))),5)];
                    sigs_mi = [sigs_mi; data((indexs(ind(lk))),4)];
                end
        end
    end
end
tot_num = numel(data(:,1));

end