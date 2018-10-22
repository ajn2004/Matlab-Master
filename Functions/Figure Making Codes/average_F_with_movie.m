%% Average F with Frames
% A script to make a movie attached with the average F obtained from user
% defined regions. The resulting movie will show a time series, location in
% time, and corresponding frame

%% User Variables
pix = 7; % pixel square radius for selecting regions

i1 = readtiff(); % load image
imagesc(std(i1,1,3)); % show standard deviation of image
num = input('How many boutons will be selected? '); % ask user for number of regions to be selected
[x,y] = ginput(num); % user selection of regions
wind = -pix:pix;
tf = [];
for i = 1:numel(x)
    sub = i1(round(y) + wind, round(x) + wind,:); % select subregion
    sf = sum(sum(sub)); % add all pixels in sub region to get flour/frame value
    tf = [tf,sf(:)];  % collimate summed flour value
end

mflour = mean(tf,2); %take the mean value over all boutons

[m,n,o] = size(i1); % get size of i1

for i = 1:o % loop over all frames
    
    % Image on the left is the raw camera image
    subplot(1,2,1);
    set_scale(i1(:,:,i),0.133,6); % scales brightness to show fluorescent response
    hold on
    draw_boxes([x,y],pix);
    title('Raw Camera image');
    axis image
    
    % Image on right is dynamic F trace
    subplot(1,2,2);
    if i < 200
        plot(t(1:300),mfluor(1:300));
        hold on
        plot([t(i),t(i)],[min(mfluor),max(mfluor)],'r');
        ylim([min(mfluor),max(mfluor)])
    else
        plot(t(i-100:i+100),mfluor(i-100:i+100))
        hold on
        plot([t(i),t(i)],[min(mfluor),max(mfluor)],'r');
        xlim([t(i-100), t(i+100)]);
        ylim([min(mfluor),max(mfluor)])
    end
    hold on
    title('average F over buotons');
    xlabel('time(s)');
    ylabel('F(au)');
    legend('F','current time');
    hold off
    drawnow;
    M(i) = getframe(gcf);
end