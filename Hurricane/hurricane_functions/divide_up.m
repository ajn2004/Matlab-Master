function [i2, fnum, cents] = divide_up(i1,pixw, dps)
[m,n,o] = size(i1);
fnum = [];
i2 = [];
cents = [];
wind = -pixw:pixw;
[X, Y] = meshgrid(wind,wind);
% [XP, YP ] = meshgrid(wind, wind);
for j = 1:o % loop over all frames
    
    [row, col] = find(dps(:,:,j) > 0); % find all peaks in Das PeakS (DPS) image
    
    
    for i = 1:numel(row) % loop over all found areas
        
        
%         if col(i) - pixw > 0 && col(i) + pixw <= n && row(i) - pixw > 0 &&  row(i) + pixw < m  % if cut out region is within the image, cut it out
            try
            % Center onto the brightest picture within a 5x5 area of the center in the reason, if you're
            % out of bounds it will catch the last error and use the
            % previous values
                [nrow, ncol] = find(i1(row(i) + wind, col(i) + wind,j) == max(max(i1(row(i) + wind, col(i) + wind,j))),1);
                col(i) = (col(i) - pixw) + ncol; % matlab's find function starts the subarray at 1 which  = col - pixw
                row(i) = (row(i) - pixw) + nrow;
                i1s = i1(row(i) + wind,col(i)+wind,j);
                xc = X.*i1s/(sum(i1s(:)));
                yc = Y.*i1s/(sum(i1s(:)));
                col(i) = col(i) - round(sum(xc(:))); % matlab's find function starts the subarray at 1 which  = col - pixw
                row(i) = row(i) - round(sum(yc(:)));

            catch lster
%                 disp(lster)
            end % end catch loop
            
            %now segment area with updated values
            if col(i) - pixw > 0 && col(i) + pixw <= n && row(i) - pixw > 0 &&  row(i) + pixw < m  % if cut out region is within the image, cut it out
                isub = i1(row(i) + wind, col(i) + wind,j); % cut out image
                i2 = cat(3, i2,isub); % linearize into array for subsequent imaging.
                fnum = [fnum;j]; % build framenum_all variable
%                 x_cm = sum(sum(isub.*XP))/(sum(isub(:)));
%                 y_cm = sum(sum(isub.*YP))/(sum(isub(:)));
                cents = [cents;[col(i), row(i)]]; % build centers found variables
            end
            
        end
 
%     fnum = [fnum;j*ones(numel(row),1)]; % build framenum_all variable
%     cents = [cents;[col, row]]; % build centers found variables
end

