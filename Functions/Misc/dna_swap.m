%Revesse compliment
% Takes in 5' 3' bottom strand, returns corresponding 5' 3' compliment
str = 'atactcgagGCCATAGCACACTTTGTACACCGggcgctagcGAAATCGGCACCGGCTTC';

for i =1:numel(str)
    n = dna_swap(str(end-i+1));
    rts(i) = n;
end
print(rts)


function o = dna_swap(i)
switch(i)
    case 'a'
        o='t';
    case 'g'
        o = 'c';
    case 'c'
        o = 'g';
    case 't'
        o = 'a';
end
end

