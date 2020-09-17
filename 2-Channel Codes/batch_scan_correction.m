current_position = pwd;
up_2_index = strfind(current_position,'\Analysis');
image_path = current_position(1:up_2_index);
files = dir('*dast_tol*');
mkdir('scan_corrected');
for i = 9
    image_file_name = [image_path, files(i).name(1:end-13),'.tif'];
    image_ruler_name = [image_path, files(i).name(1:end-14),'scan.tif'];
    file_list = {files(i).name,image_file_name,image_ruler_name};
    t(i) = laser_scan_correction(file_list);
    disp(['Completed ', num2str(i), ' files in ', num2str(sum(t)),'s'])
    disp(['Estimating completion in ', num2str(mean(t)*(numel(files)-1)),'s'])
end
