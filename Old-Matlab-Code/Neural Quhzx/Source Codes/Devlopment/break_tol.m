% % A function to apply weak tolerances to breakface localizations
% function break_tols(fname) 
% load(fname);
clearvars
load('hello.mat');
varys = {'xf_all','yf_all','N','sigx_all','sigy_all', 'off_all','off_crlb','xf_crlb','yf_crlb','N_crlb','sigy_crlb','sigx_crlb','llv'};

% Remove all negative crlbs (they are variances and should be inherently
% positive
ind = find(N_crlb <= 0 | off_crlb <= 0 | sigx_crlb <= 0 | sigy_crlb < 0 | xf_crlb <= 0 | yf_crlb <= 0);

for i = 1:numel(varys)
    eval([varys{i},'(ind) = [];']);
end


% Build Fractional Uncertanties
frn = N_crlb.^0.5./N;
fro = off_crlb.^0.5./off_all;
frx = sigx_crlb.^0.5./sigx_all;
fry = sigy_crlb.^0.5./sigy_all;
varys = {'frn','fro','frx','fry','xf_all','yf_all','N','sigx_all','sigy_all', 'off_all','off_crlb','xf_crlb','yf_crlb','N_crlb','sigy_crlb','sigx_crlb','llv'};

bigf = [llv, frn, fro, frx, fry, xf_crlb, yf_crlb];
bigf = bigf - mean(bigf);
vof =cov(bigf);
[D,I] = eig(vof);
sdf = D(:,end).'*bigf.';

% bigv = [N, sigx_all, sigy_all, off_all, xf_all, yf_all];
% bigv = bigv - mean(bigv);
% vov =cov(bigv);
% [Dv,I] = eig(vov);
% sdv = Dv(:,end).'*bigv.';

% bigc = [N_crlb, sigx_crlb, sigy_crlb, off_crlb, xf_crlb, yf_crlb];
% bigc = bigc - mean(bigc);
% voc =cov(bigc);
% [Dc,I] = eig(voc);
% sdc = Dc(:,end).'*bigc.';

ind1 = sdf < 0;
ind2 = sdf > 0;
subplot(1,2,1); plot(xf_all(ind1),yf_all(ind1),'.b'); title('Less than 0');
subplot(1,2,2); plot(xf_all(ind2),yf_all(ind2),'.b'); title('Greater than 0');
while true
    inp = input('(H)igh or (L)ow?','s');
    if strcmp(inp,'l') || strcmp(inp,'L')
        ind = ind2;
        break
    elseif strcmp(inp,'h') || strcmp(inp,'H')
        ind = ind1;
        break
    else
    end
end
close all
for i = 1:numel(varys)
    eval([varys{i},'(ind) = [];']);
end

ind = off_all < 0 | sigx_all < 1.1 | sigx_all > 6 | sigy_all < 1.1 | sigy_all > 6 | sigy_crlb > 10 | sigx_crlb >10 | xf_crlb > 1 | yf_crlb > 1 | frn >= 1 | fro >=1 | frx >=1 | fry >= 1;


for i = 1:numel(varys)
    eval([varys{i},'(ind) = [];']);
end
% plot(xf_all,yf_all,'.b')
histogram((sigx_all.*sigy_all).^0.5);
% plot(sigx_all, sigy_all,'.b')