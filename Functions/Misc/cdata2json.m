function output_string = cdata2json(cdata)
% We're changing the x-y-z cdata to a json string

% data curation
% transform orange onto x
% [xfo, yfo] = make_nn_channel_transform(cdata.orange.xf,cdata.orange.yf);
xfo = cdata.orange.xf;
yfo = cdata.orange.yf;
% Readjust image so that mean of red distribution is (0,0, arbitrary)
xfo = xfo - mean(cdata.red.xf);
yfo = yfo - mean(cdata.red.yf);
zfo = cdata.orange.zf - mean(cdata.red.zf);
xfr = cdata.red.xf - mean(cdata.red.xf);
yfr = cdata.red.yf - mean(cdata.red.yf);
zfr = cdata.red.zf - mean(cdata.red.zf);


output_string = '{ "red": { "x": [';
for i = 1:numel(cdata.red.xf)
    output_string = [output_string , num2str(xfr(i)), ','];
end
output_string = [output_string(1:end-2), '], "y": ['];
for i = 1:numel(cdata.red.xf)
    output_string = [output_string , num2str(yfr(i)), ','];
end
output_string = [output_string(1:end-2), '], "z": ['];
for i = 1:numel(cdata.red.xf)
    try
%         output_string = [output_string , num2str(cdata.red.zf_raw(i) - mean(cdata.red.zf_raw)), ','];
    output_string = [output_string , num2str(zfr(i)), ','];
    catch
        output_string = [output_string , num2str(cdata.red.zf(i)), ','];
    end
end
output_string = [output_string(1:end-2), ']}, "orange": { "x" : ['];
for i = 1:numel(cdata.orange.xf)
    output_string = [output_string , num2str(xfo(i)), ','];
end
output_string = [output_string(1:end-2), '], "y": ['];
for i = 1:numel(cdata.orange.xf)
    output_string = [output_string , num2str(yfo(i)), ','];
end
output_string = [output_string(1:end-2), '], "z": ['];
for i = 1:numel(cdata.orange.xf)
    try
%         output_string = [output_string , num2str(cdata.orange.zf_raw(i)- mean(cdata.red.zf_raw)), ','];
output_string = [output_string , num2str(zfo(i)), ','];
    catch
        output_string = [output_string , num2str(cdata.orange.zf(i)), ','];
    end
end
output_string = [output_string(1:end-2), ']}}'];