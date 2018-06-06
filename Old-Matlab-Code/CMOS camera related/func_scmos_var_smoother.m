function [smoothed_im] = func_scmos_var_smoother(i1,sig,varian)

index1 = logical(zeros(numel(i1(:,1)),numel(i1(1,:))));
index2 = index1;
S1 = zeros(numel(i1(:,1)),numel(i1(1,:)));
S2 = S1;
for i = floor(2*sig+1):numel(i1(:,1))-floor(2*sig+1)
    for j = floor(2*sig+1):numel(i1(1,:))-floor(2*sig+1)
        index1(i-sig:i+sig,j-sig:j+sig) = 1;
        index2(i-2*sig:i+2*sig,j-2*sig:j+2*sig) = 1;
        counts1 = double(i1(index1));
        counts2 = double(i1(index2));
        vari1 = double(varian(index1));
        vari2 = double(varian(index2));
        S1(i,j) = sum(counts1./vari1)/sum(1/vari1);
        S2(i,j) = sum(counts2./vari2)/sum(1/vari2);
        index1 = logical(zeros(numel(i1(:,1)),numel(i1(1,:))));
        index2 = index1;
        imshow(S1)
    end
%     waitforbuttonpress;
end
smoothed_im = S1-S2;
smoothed_im = smoothed_im.*(smoothed_im>0);