function cdata = model_drift_correction(cdata, align_color, frames_in_chunk )

% Build initial coordinate set and overlay via neural net transform
% xfo and xfr were initially to indicate orange from red, but given the
% current analysis is color indiscriminante this is a vestige from
% development that will be kept as xfr: reference distribution xfo: off
% color distribution
[xfo, yfo] = make_nn_channel_transform(cdata.orange.xf,cdata.orange.yf);

if strcmpi(align_color,'orange') % If we specify orange, make reference distribution be orange molecules
    [xfr, yfr] = make_nn_channel_transform(cdata.orange.xf,cdata.orange.yf);
    zfr = -cdata.orange.zf;
    xfo = cdata.red.xf;
    yfo = cdata.red.yf;
    zfo = -cdata.red.zf;
    frame_o = cdata.red.framenumber;
    frame_r = cdata.orange.framenumber;
    xfo = xfo-mean(xfr);
    zfo = zfo-mean(zfr);
    yfo = yfo-mean(yfr);
    xfr = xfr-mean(xfr);
    zfr = zfr-mean(zfr);
    yfr = yfr-mean(yfr);
else
    [xfo, yfo] = make_nn_channel_transform(cdata.orange.xf,cdata.orange.yf);
    zfo = -cdata.orange.zf;
    xfr = cdata.red.xf;
    yfr = cdata.red.yf;
    zfr = -cdata.red.zf;
    xfo = xfo-mean(xfr);
    zfo = zfo-mean(zfr);
    yfo = yfo-mean(yfr);
    xfr = xfr-mean(xfr);
    zfr = zfr-mean(zfr);
    yfr = yfr-mean(yfr);
    frame_r = cdata.red.framenumber;
    frame_o = cdata.orange.framenumber;
end

% Make data arrays
r_data = [xfr, yfr, zfr];
o_data = [xfo, yfo, zfo];

% Build t0 model by segmenting out early data
index = frame_r >= 0 & frame_r < frames_in_chunk;
data0 = r_data(index,:);