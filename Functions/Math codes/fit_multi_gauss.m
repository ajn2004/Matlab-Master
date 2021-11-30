function [beta, cost] = fit_multi_gauss(data, N)

[N_red,edges] = histcounts(data,round(numel(data)^(2/3)));
for i  = 1:numel(edges)-1
    cent_red(i) = (edges(i) + edges(i+1))/2;
end

if N < 1
    N = 1;
end

% Generate N-Gaussian function w/ guess parameters
string = 'modelfun = @(b,x)(';
beta_guess = 'beta0 = [max(N_red), sum(cent_red.*N_red./sum(N_red)), 0.02, ';
for i = 1:N
    term = ['b(',num2str(1 + (i-1)*3),')*exp(-((x-b(', num2str(2 + (i-1)*3),')).^2)/(2*b(',num2str(3 + (i-1)*3),')^2))'];
    if i > 1
        string = [string,' + ',term];
        beta_guess = [beta_guess, num2str(1/i), '*max(N_red), sum(cent_red.*N_red./sum(N_red)) - ',num2str(0.01*(i-1)),', ',num2str(i),'*0.02,'];
    else
        string = [string,term];
    end
end

% Finishing Touch on strings
string = [string,' + b(', num2str(i*3+1),'));'];
beta_guess = [beta_guess,'0]'];

% evaluate strings
eval(string);
eval(beta_guess);

% perform fit
beta = nlinfit(cent_red,N_red,modelfun,beta0);
% Extract Cost
cost = sum(abs(modelfun(beta,cent_red) - N_red));
