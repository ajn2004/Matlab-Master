function makebox(cents,psize)
plot([cents(1) - psize, cents(1) + psize] , [cents(2) - psize, cents(2) - psize], 'c')
plot([cents(1) + psize, cents(1) + psize] , [cents(2) - psize, cents(2) + psize], 'c')
plot([cents(1) - psize, cents(1) + psize] , [cents(2) + psize, cents(2) + psize], 'c')
plot([cents(1) - psize, cents(1) - psize] , [cents(2) + psize, cents(2) - psize], 'c')
end