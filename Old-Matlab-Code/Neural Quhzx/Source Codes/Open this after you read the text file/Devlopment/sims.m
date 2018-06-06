clearvars
[xg, yg] = meshgrid(-3:3,-3:3);

Na = 500;
sigt = 1.75;
count = 1;
for off = 0:0.01:8
mols = 1000;
i1 = [];
for i = 1: mols
    z = Na/(2*pi*sigt.^2)*exp(-((xg - randn).^2 + (yg - randn).^2)/(2*sigt^2)) + off;
    %     z = 500/(2*pi*1.47*1.7)*exp(-xg.^2/(2*1.47^2)).*exp(- yg.^2/(2*1.7^2));
    zi = single(imnoise(uint16(z), 'Poisson'));
    i1 = [i1,zi(:)];
end

[xf_all,xf_crlb, yf_all,yf_crlb, N,  N_crlb,  sigx,  sigx_crlb,sigy,    sigy_crlb,off_all, off_crlb, llv] = full_chain_loc(single(i1),100);
t(count) = sum(llv<0);
offs(count) = off;
count = count+1;
end
plot(offs, t/mols)