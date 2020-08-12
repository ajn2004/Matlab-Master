function mW = power_converter(mA,lamb)

% Emperical data for 
pows = [0,5,10,25,50,75,100];
f88 = [0, 0.817 1.621 4.119 7.412 10.224 12.686];
f61 = [0, 1.656 2.628 6.269 10.959 14.927 18.260];
s37 =[0, 1.483 2.919 7.240 12.478 16.720 21.087];
mW = mA*0;
m = numel(mA);
if m ~= numel(lamb)
    lamb = ones(m,1)*lamb(1);
end
for i = 1:m
switch lamb(i)
    case 488
        mW(i) = spline(f88,pows,mA(i));
    case 561
        mW(i) = spline(f61,pows,mA(i));
    case 637
        mW(i) = spline(s37,pows,mA(i));
    otherwise
        mW(i) = -1;
end
end