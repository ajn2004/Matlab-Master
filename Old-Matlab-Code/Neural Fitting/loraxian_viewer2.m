%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Loraxian viewer
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i1 = zeros(9);
[X, Y] = meshgrid(-4:4,-4:4);
[fname, fpath] = uigetfile('*.mat');
load([fpath,fname]);
figure('Units','Normalized','Outerposition', [0 0 1 1]);
histnn =[];
histcc =[];
histNn= [];
histBn =[];
histsx =[];
histNm = histsx;
histBm = histNm;
histsxm = histsx;
for i =1:10000
    x0 =0.5*randn;
    y0 = 0.5*randn;
    N = randi(800)+50;
    sigx = 2*rand+0.5;
    sigy = sigx;
    B = rand*(6);
    
    % make it gauss
    i1x = 1/2.*(erf((X - x0 + 1/2)./(2*sigx^2)^0.5)-erf((X - x0 - 1/2)./(2*sigx^2)^0.5)); % error function of x integration over psf for each pixel
    i1y = 1/2.*(erf((Y - y0 + 1/2)./(2*sigy^2)^0.5)-erf((Y - y0 - 1/2)./(2*sigy^2)^0.5)); % error function of y integration over psf for each pixel
    i1 =  N * i1x.*i1y+B;
    i2 = double(imnoise(uint16(i1), 'poisson'))+.00001;
    % Neural Calculation
    X1 = [ 1, i2(:).'];
    z2 = X1*theta1.';
    a2 = [ones(numel(z2(:,1)),1), sigmoid(z2)];
    z3 = a2*theta2.';
    a3 = sigmoid(z3);
    
    %unscale
    outys = a3.*ysc + ymi;
    xf = outys(1);
    yf = outys(2);
    Nf = outys(3);
    sxf = outys(4);
    syf = outys(5);
    Bf = outys(6);
    
    xc =0;
    yc =0;
    for j = 1:9
        yc = yc + j*sum(i2(j,:));
    end
    for k = 1:9
        xc = xc + k*sum(i2(:,k));
    end
    yc = yc / sum(i2(:));
    xc = xc / sum(i2(:));
    [xm, ym,Nm, sigxm, sigym, offset] = func_mle(i2, xc-5, yc-5);
    histNn = [histNn;Nf-N];
    histNm = [histNm;Nm-N];
    histBn = [histBn;Bf-B];
    histBm = [histBm;offset-B];
    histsx = [histsx;sxf - sigx];
    histsxm = [histsxm;sigxm - sigx];
    %represent
    subplot(3,2,[1,3]);
    imagesc(i2); colormap('gray');
    hold on
    plot(xf+5,yf+5,'.b','MarkerSize',20);
    
    plot(x0+5,y0+5,'.r','MarkerSize',20);
    plot(xm+5,ym+5,'.g', 'MarkerSize',20);
    hold off
    subplot(3,2,2); histogram(histNn,-150:5:150);hold on; histogram(histNm,[-150:5:150]);hold off
    subplot(3,2,4); histogram(histBn,-5:0.1:5);hold on; histogram(histBm,[-5:0.1:5]);hold off
    
    
    drawnow
    histnn = [histnn;((xf-x0)^2 +(yf-y0)^2)^0.5];
    histcc = [histcc;((xm-x0)^2 +(ym-y0)^2)^0.5];
    subplot(3,2,5); histogram(histnn,[0:0.005:.75]);
    hold on
    histogram(histcc,[0:0.005:0.75]);
    hold off
    legend('NN','MLE')
    title('Histogram of diffrences')
    subplot(3,2,6); histogram(histsx,-0.4:0.005:0.4); hold on;histogram(histsxm,[-0.4:0.005:0.4]); title('Histogram of difference in sigx'); hold off
    pause(0.5)
end