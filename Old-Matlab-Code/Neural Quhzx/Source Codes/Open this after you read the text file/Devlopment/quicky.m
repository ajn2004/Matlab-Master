[xg, yg] = meshgrid(-3:3,-3:3);
z = 500/(2*pi*1.75*1.75)*exp(-(xg.^2 + yg.^2)/(2*1.75^2));
imagesc(z)
[xf_all,xf_crlb, yf_all,yf_crlb, N,  N_crlb,  sigx,  sigx_crlb,sigy,    sigy_crlb,off_all, off_crlb, llv] = full_chain_loc([single(z(:)*0),single(z(:))],10);