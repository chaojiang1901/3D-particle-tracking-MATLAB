function [inTracksStruct,outTracksStruct] = SelectTracksInRegion(TracksStruct,Region)

for itrack=1:length(TracksStruct)
points = TracksStruct(itrack,1).points;
mean_x(itrack,1) = mean(points(:,1));
mean_y(itrack,1) = mean(points(:,2));
end
ROI_xv=Region(:,1);ROI_yv=Region(:,2);
in = inpolygon(mean_x,mean_y,ROI_xv,ROI_yv);
inTracksStruct = TracksStruct(in);
outTracksStruct=TracksStruct(~in);
