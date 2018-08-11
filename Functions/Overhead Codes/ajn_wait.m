function ajn_wait(t, i, o)
clc;
% here t is the time it takes to complete 1 iteration, I is the ith
% iteration and o is the total number of iterations needed to complete
    if numel(t) > 1 % only begin analysis if there is more than 1 t value
        t1(1) = t(1);
        for j = 2:numel(t)
            t1(j) = t(j) - t(j-1); % build a vector of differences
        end
    tm_rem = mean(t)*(o-i);  % average time/element value multiplied remaining elements
    tm1_rem = sum(t)/i*(o-i); % total time / all elements multiplied by remaining elements
    tm_rem = (tm_rem + tm1_rem)/2;
    sec_rem = round(mod(tm_rem,60));
    mins = floor(tm_rem/60);
    hr_rem = floor(mins/60);
    min_rem = mod(mins,60);
%     clc
    disp([num2str(100*i/o),'% complete']);
    disp(['Time remaining is ~', num2str(hr_rem),' hrs, ', num2str(min_rem),' mins, ', num2str(sec_rem),'s']);
    end