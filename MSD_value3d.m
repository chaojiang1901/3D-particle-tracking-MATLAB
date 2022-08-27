function [ msd ] = MSD_value3d(window_frameNum,window_points, msd_size)
windowlength=length(window_frameNum); 
temp_frameNum= repmat(window_frameNum,1,windowlength);
temp1_frameNum= temp_frameNum';
temp_x = repmat(window_points(:,1),1,windowlength);
temp1_x = temp_x';
temp_y = repmat(window_points(:,2),1,windowlength);
temp1_y = temp_y';%×ªÖÃ
temp_z = repmat(window_points(:,3),1,windowlength);
temp1_z = temp_z';%×ªÖÃ

frame_interval = temp_frameNum-temp1_frameNum;
frame_displacement = ((temp_x-temp1_x).^2+(temp_y-temp1_y).^2+(temp_z-temp1_z).^2);

for imsd=1:msd_size
    msd(imsd,1) = mean(frame_displacement(frame_interval==imsd));
end


