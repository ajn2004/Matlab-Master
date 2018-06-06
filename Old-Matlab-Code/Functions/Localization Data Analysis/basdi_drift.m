% function [xf_dc,yf_dc, s] = basdi_drift(xf_all,yf_all,framenum_all, q)
% This function will allow me to convert our localization data to the BaSDI
% approach described in 2015 BJ, basdi codes are from authors
% bx = ceil(max(xf_all));
% by = ceil(max(yf_all));
% mx = floor(min(xf_all));
% my = floor(min(yf_all));

w = round((bx + 10)*q*1000/30);
h = round((by + 10)*q*1000/30);
% xf = xf_all;
% yf = yf_all;
clear O
for i = 1:max(framenum_all)
    clear o
    id = find(framenum_all == i);
    if ~isempty(id)
%         disp(id)
    o(:,2) = round(xf_all(id)*q*1000/30);
    o(:,1) = round(yf_all(id)*q*1000/30);
    O{i} = o;
    else
        O{i} = [];
    end
    
end

[s] = BaSDI_main(O,h,w);
xf_dc = [];
yf_dc = [];
% for i = 1:numel(d(:,1))
%     id = find(framenum_all == i);
%     xf_dc = [xf_dc;xf_all(id) - d(i,1)];
%     yf_dc = [yf_dc;yf_all(id) - d(i,2)];
% end