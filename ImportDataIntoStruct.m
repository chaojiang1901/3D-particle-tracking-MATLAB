% read data.txt files and build TracksStruct

 function [Struct] = ImportDataIntoStruct(TrackLength,x)

[FileName,PathName] = uigetfile('*.txt','Select the data file','MultiSelect','on');

if iscell(FileName)
for ifile = 1:length(FileName)
file{ifile}=strcat(PathName,FileName{ifile});
fid=fopen(file{ifile});
itrack=1;
judge=1;
clear TracksStruct
TracksStruct(1,1) = struct('frameNum',[],'points',[],'intensities',[]);
while judge
    target=fread(fid,2,'*char');
    target=target';
    judge=strcmp(target,'%%');
       if judge
           fgetl(fid);
           trackdata=textscan(fid,'%f %f %f %f %f %s');
     if length(trackdata{1})>=TrackLength && x>=length(trackdata{1});
           TracksStruct(itrack,1).frameNum = trackdata{1};
           TracksStruct(itrack,1).points = [trackdata{3},trackdata{2}]; %change the xy position
           TracksStruct(itrack,1).intensities = [trackdata{4},trackdata{5}];         
           itrack=itrack+1;                   
     end           
       end          
end         
fclose(fid);
Struct(ifile,1).TracksStruct=TracksStruct;
Struct(ifile,1).FileName = FileName{ifile};
end

else
file=strcat(PathName,FileName);
fid=fopen(file);
itrack=1;
TracksStruct(1,1) = struct('frameNum',[],'points',[],'intensities',[]);
judge=1;
while judge
    target=fread(fid,2,'*char');
    target=target';
    judge=strcmp(target,'%%');
       if judge
           fgetl(fid);
           trackdata=textscan(fid,'%f %f %f %f %f %s');
     if length(trackdata{1})>=TrackLength && x>=length(trackdata{1});
           TracksStruct(itrack,1).frameNum = trackdata{1};
           TracksStruct(itrack,1).points = [trackdata{3},trackdata{2}]; %change the xy position
           TracksStruct(itrack,1).intensities = [trackdata{4},trackdata{5}];         
           itrack=itrack+1;                   
     end           
       end          
end         
fclose(fid);
Struct(1,1).TracksStruct = TracksStruct;
Struct(1,1).FileName = FileName;
end
