multips = [];
primes = [];
for i = 2:100000
    if sum(multips == i) == 0
        primes = [primes;i];
        multips = [multips; (1:round(100000/i)).'*i];
    else
    end
end