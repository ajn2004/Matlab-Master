function index = localization_crlb_filtering(cdata,q)
% This function will filter localizations based on Carmer-Rao Lower bounds

color = 'red';
sigx_o = cdata.(color).crlbs(:,4).^0.5*q;
sigy_o = cdata.(color).crlbs(:,5).^0.5*q;
srxy = (sigx_o.^2 + sigy_o.^2).^0.5; % Quadratic combination of widths
index.red = srxy<0.015; % This has proven to be a reasonable constraint so we can hard code it into the function as opposed to making it a variable
color = 'orange';
sigx_o = cdata.(color).crlbs(:,4).^0.5*q;
sigy_o = cdata.(color).crlbs(:,5).^0.5*q;
soxy = (sigx_o.^2 + sigy_o.^2).^0.5;
index.orange = soxy<0.015;

end