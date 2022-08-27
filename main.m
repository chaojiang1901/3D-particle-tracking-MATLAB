%% import txt data from Particle tracker to matlab;
clear
ViscoCellData = ImportDataIntoStruct(50,6000); %'tracklength for filter'
%% %remove the tracks out of the region of interest;
for iSample = 1:length(ViscoCellData)
    iSample
    TracksStruct = ViscoCellData(iSample,1).TracksStruct;  
    temTracksStruct=TracksStruct;
    [inTracksStruct,outTracksStruct] = SelectTracksInRegion(TracksStruct,Region); %针对points1进行筛选
    %ViscoCellData(iSample,1).TracksStruct=inTracksStruct;
    ViscoCellData(iSample,1).TracksStruct=inTracksStruct;
end

%% %enter the coordinates of corresponding points to calculate the transmatrix;
pk1=[197.0000	49.0000
50.0000	163.0000
38.0000	330.0000];
pk2=[184.0000	44.0000
40.0000	163.0000
33.0000	333.0000];
right=pk1;
left=pk2;
right=[right,ones(3,1)]'; % right is from focused image，left is from defocused image;
left=[left,ones(3,1)]'; 
trans=left*pinv(right) %trans is the transmatrix between focused and defocused image;
%% %Convert the xy coordinates of the focus plane to the one of the defocus plane;
for iSample = 1:length(ViscoCellData)
    iSample;
    TracksStruct= ViscoCellData(iSample,1).TracksStruct; 
    for itrack =1:length(TracksStruct)
        points = TracksStruct(itrack,1).points;
        pointslength=length(points);
        x3=ones(pointslength,1);
        points3=[points,x3];
        TracksStruct(itrack,1).points3=points3;
    end
    ViscoCellData(iSample,1).TracksStruct=TracksStruct;       
end

for iSample = 1:length(ViscoCellData)
    iSample;
    TracksStruct= ViscoCellData(iSample,1).TracksStruct;     
    for itrack =1:length(TracksStruct)
        points4=[];
        points3 = TracksStruct(itrack,1).points3;
        pointslength=length(points3);
        %pointslength=2;
        for i=1:pointslength            
            xj=points3(i,:,:)';
            yj=trans*xj;
            zj=yj';
            points4=[points4;zj];
        end
        points4(:,3:m)=[];
        %x3=ones(pointslength,1);
        %points3=[points,x3];
        TracksStruct(itrack,1).points4=points4;
    end
    ViscoCellData(iSample,1).TracksStruct=TracksStruct;   
end
2
%% 3D 可以针对多条轨迹进行集中提取图片
str='K:\3d-tracking-Data-2020\20220121-3D-Ttracking-isotropic-T40dextran\JC17\def\';
path='K:\3d-tracking-Data-2020\20220121-3D-Ttracking-isotropic-T40dextran\JC17\image\';
%path2='G:\experiment\20170120-dual-3d\fish11D-defocus2-nouseing-jc';
mi=37;%define the size of ROI of diffraction ring;
r=(mi+1)/2;
m2=30;
%r1=round(r*1.5);
x=1:mi;
y=1:mi;
[X,Y]=meshgrid(x,y);
[TH,R]=cart2pol(X-r,Y-r);
for iSample = 1:length(ViscoCellData)
    iSample;
    TracksStruct= ViscoCellData(iSample,1).TracksStruct; 
    for itrack =1:length(TracksStruct);
    %for itrack=3:5
        points4 = TracksStruct(itrack,1).points4;
        frameNum = TracksStruct(itrack,1).frameNum;
        ww=0;
        for i=1:length(frameNum);
            ww=ww+1;
            %str='G:\experiment\20170120-dual-3d\201701120-dual-fish11cut\201701120-dual-fish11cut';
            I=imread([str,num2str(frameNum(i),'%04u'),'.tif']);
            %I=im2double(I);
            %imagebp=bpass(image,0.1,13);
            %pk=pkfnd(imagebp,0.001,13);
            [M,N]=size(I);
            I1=zeros(M+2*m2,N+2*m2);
            I1(m2+1:M+m2,m2+1:N+m2)=I;
            %j=ceil(points4(i,1)+1);
            j=round(points4(i,1)+1);
            %k=ceil(points4(i,2)+1);
            k=round(points4(i,2)+1);
            j1=j+m2;k1=k+m2;
            try          
            Z=I1(k1-r+1:k1+r-1,j1-r+1:j1+r-1);
            Z=mat2gray(Z);
                img=Z;
                filename=strcat('traj',num2str(itrack,'%04u'),'--frame',num2str(frameNum(i),'%04u'),'.tif');
                %filename=strcat(
                pathfile=fullfile(path,filename);
                imwrite(img,pathfile,'tif');
            %else
            %    img=Z;
            %    filename=strcat(num2str(frameNum(i)),'.tif')
            %    pathfile=fullfile(path2,filename);
            %    imwrite(img,pathfile,'tif');
            %    disp(strcat('Brightcontrast is bad: ',num2str(frameNum(i))))
            catch
            end
            %end    
        end
    end
end
%%  calculate the radius of ring
% 
%clear aResultStruct;
%frameNum = TracksStruct.frameNum;
%str2='K:\3d-tracking-Data-2020\20220117-3D-Ttracking-isotropic-T40dextran-Glucose\JC15-glucose\image\';
str2=path;
%str3='D:\20171009\7\traj5962\';
%mi=37; % 
rmi=(mi+1)/2;
%rmi1=rmi+5; %
x=1:mi;
y=1:mi;
[X,Y]=meshgrid(x,y);
[X1,Y1]=meshgrid(-rmi+1:1:rmi-1);
[TH,Rall]=cart2pol(X-rmi,Y-rmi);
IJ=unique(Rall);
for iSample = 1:length(ViscoCellData)
    iSample;
    TracksStruct= ViscoCellData(iSample,1).TracksStruct; 
    for itrack=1:length(TracksStruct)
    %for itrack=56
        points = TracksStruct(itrack,1).points; %用于后续计算，在这个提取半径的部分，不必要。
        frameNum = TracksStruct(itrack,1).frameNum;
        intensities=TracksStruct(itrack,1).intensities;
        if  ~isempty(frameNum)   
        %points=points4;
        ix=1;
        iy=length(frameNum);     
        %iy=round(length(frameNum)/2); 
        %iy=9;
        aRlast=[];
        for i=ix:iy;
            filename2=strcat('traj',num2str(itrack,'%04u'),'--frame',num2str(frameNum(i),'%04u'),'.tif');
            %filename=strcat(num2str(i-1,'%04u'),'.tif');
            %filename2=strcat(num2str(frameNum(i),'%04u'),'.tif');
            oanvn=exist([str2,filename2]);
            if oanvn>0
                    %IU=imread([str3,num2str(i-1,'%04u'),'.tif']); %从零开始，减去初始帧数。
                    IU=imread([str2,filename2]);
                    IU=mat2gray(IU);
                    IU4=wiener2(IU,[3,3]);
                    IU1=bpass(IU4,0.01,3);
                    %IU2=IU1;
                    IU2=3*IU1;
                    BW=im2bw(IU2,2*mean2(IU2));
                    rmin=1;
                    rmax=17;
                    P=4; % precision (large P will use a huge amount of memory)
                    FS=4; % Gaussian filter size (large FS will be extremely slow)
                    [xc,yc,rh]=circleHough(BW,rmin,rmax,P,FS); %https://ww2.mathworks.cn/matlabcentral/fileexchange/46380-circle-hough-transform?s_tid=FX_rc1_behav
                    %i;
                    XXrsquare=[];
                    YYsse=[];
                    Rrrr=[];
                %if rh<10
                    %for b=rh;    
                    for b=rh-3:3:rh+3
                        try
                            [fitresult, gof] = createFitradius(Y1,X1,IU1,b); 
                            XXrsquare=[XXrsquare,gof.rsquare];
                            YYsse=[YYsse,gof.sse];
                            Rrrr=[Rrrr,fitresult.p]; 
                        catch                        
                            disp(strcat('bad '));
                        end                        
                    end     
                    if  ~isempty(XXrsquare)                    
                    [idxrsquare,idxrsquareloc]=max(XXrsquare);
                    Rfinal=Rrrr(idxrsquareloc);
                    rsquarefinal=idxrsquare;
                    ssefinal=YYsse(idxrsquareloc);
                    if idxrsquare>0.8 %& Rfinal<rh+2 & rh-2<Rfinal
                        %if rh>6.2;
                        aRlastStruct=[frameNum(i),points(i,1),points(i,2),Rfinal,intensities(i,1),rh,rsquarefinal,ssefinal];
                        %else 
                        %aRlastStruct=[frameNum(i),points(i,1),points(i,2),Rfinal,intensities(i,1),rh,rsquarefinal,ssefinal];
                        %end
                %clear rowsHer; clear colsHer;
                        aRlast=[aRlast;aRlastStruct];
                    else
                        %disp(strcat('fitting is bad --',num2str(i)));
                    end     
                    else
                        %disp(strcat('This image is bad--',num2str(i)));
                    end
                %else
                %    disp(strcat('pass this image,the noise is so high ','-',num2str(i,'%04u')))            
                %end
            else
                %disp(strcat('pass this image,the noise is so high ','-',num2str((itrack),'%04u'),'-',num2str(frameNum(i),'%04u')))
                %disp(strcat('pass this image,the noise is so high ','-',num2str(i,'%04u')));
            end 
        
        end
        else 
            aRlast=[];
        end
        TracksStruct(itrack,1).aRlast=aRlast;
        clear aResultStruct;
    end
    
end
ViscoCellData(iSample,1).TracksStruct=TracksStruct;
%% %select tracks longer than 50 frames in a new TracksStruct;
clear TracksStructnew;
for iSample = 1:length(ViscoCellData)
    iSample
    TracksStruct = ViscoCellData(iSample,1).TracksStruct;  
    n=7;
    m=0;
    for itrack =1:length(TracksStruct)
    %for itrack=n:n
        aRlast = TracksStruct(itrack,1).aRlastnew; 
        %aRlastnew=aRlast(aRlast(:,4)>Zmiddle-Zheight & aRlast(:,4)<Zmiddle+Zheight,:); 
            if length(aRlast)>50               %aRlastnew=aRlast(aRlast(:,4)<10,:);
               m=m+1;
                TracksStructnew(m,1).aRlastnew=aRlast;  
                TracksStructnew(m,1).Tracknum=itrack;  
            else
            end
        end
    end
    ViscoCellData(iSample,1).TracksStructnew=TracksStructnew; 
%% Get the 3D coordinates (X,Y,Z)
TimeResolution = 0.03;
MSD_length = 33;
limit_length = 50;
kk=0.235; % um the calibration between Z coordinates and the ring radius;
iSample=1;
nSample = length(ViscoCellData);
for iSample = 1:nSample
TracksStructnew = ViscoCellData(iSample,1).TracksStructnew;  
MSD_select=[];
MSD_select1=[];
iMSD=0;
clear Track;
for itrack = 1:length(TracksStructnew)
    aRlastnew= TracksStructnew(itrack,1).aRlastnew;
    frameNum = aRlastnew(:,1);
    points = aRlastnew(:,2:4);
    intensities=aRlastnew(:,5);
    points(:,1)=points(:,1)*0.2667;
    points(:,2)=points(:,2)*0.2667;
    points(:,3)=points(:,3)*kk;
    Track(itrack,1).frameNum=frameNum;
    Track(itrack,1).points=points;
    Track(itrack,1).intensities=intensities;
end
ViscoCellData(iSample,1).Track=Track;
clear Track;
end
1
%% Calculate MSD and fitting D and alpha from several cells;
%y=d*x^alpha
% jc
limit_length=50;MSD_length=33;
TimeResolution = 0.03;
timelag=TimeResolution;
alphasize=33;
MSDnew=[];
MSDnew1=[];
MSDnew2=[];
MSDnewx=[];
MSDnewy=[];
MSDnewz=[];
for iSample=1:length(ViscoCellData1)
    ViscoCellData=ViscoCellData1(iSample,1).ViscoCellData;
    %Track=ViscoCellData.Track;
    Track=ViscoCellData.Track;
%      Track=ViscoCellData.DirectStruct;
%    Track=ViscoCellData.DirectwholeStruct;    
%     Track=ViscoCellData.DiffStruct;
%     Track=ViscoCellData.DiffusionStruct;
%        Track=ViscoCellData.SubDiffStruct;
%        Track=ViscoCellData.ConfinedStruct;
     MSD_select=[];
        MSD_select1=[];
        MSD_select2=[];
        MSD_selectx=[];
        MSD_selecty=[];
        MSD_selectz=[];
    MSD_measureTime = [];
    iMSD=0;
    %for itrack=8;
    if ~isempty(Track(1).frameNum);
    for itrack = 1:length(Track)
       
        frameNum=Track(itrack,1).frameNum;
        points=Track(itrack,1).points;
        %intensities=Track(itrack,1).intensities;
        if length(points(:,1))>limit_length
           iMSD=iMSD+1;
         MSD_select(:,iMSD) = MSD_value(frameNum,points(:,1:2), MSD_length);%*0.2667*0.2667;
         MSD_select1(:,iMSD) = MSD_value3d(frameNum,points(:,1:3), MSD_length);%*0.2667*0.2667*0.18;
         MSD_select2(:,iMSD) = MSD_value1d(frameNum,points(:,3), MSD_length);
         MSD_selectx(:,iMSD) = MSD_value1d(frameNum,points(:,1), MSD_length);%*0.2667*0.2667;
         MSD_selecty(:,iMSD) = MSD_value1d(frameNum,points(:,2), MSD_length);%*0.2667*0.2667*0.18;
         MSD_selectz(:,iMSD) = MSD_value1d(frameNum,points(:,3), MSD_length);
        end
    end
    MSDnew=[MSDnew,MSD_select]; %0-xy-2D
    MSDnew1=[MSDnew1,MSD_select1];%1-xyz-3D
    MSDnew2=[MSDnew2,MSD_select2];%2-z
    MSDnewx=[MSDnewx,MSD_selectx];
    MSDnewy=[MSDnewy,MSD_selecty];
    MSDnewz=[MSDnewz,MSD_selectz];

98889
% %% %每条轨迹计算Dx Dy Dz Dxyz Dxy alpha
%for iSample = 1:length(ViscoCellData1)
%ViscoCellData= ViscoCellData1(iSample,1).ViscoCellData; 
%cellcenter= ViscoCellData1(iSample,1).cellcenter; 
    %Track=ViscoCellData.Track;
%     DaF4L=[];Da1=[];
%     for itrack = 1:length(Track)
%             frameNum=Track(itrack,1).frameNum;
%             points=Track(itrack,1).points;
%             %intensities=Track(itrack,1).intensities;
%             if length(points(:,1))>limit_length
%             meanx=mean(points(:,1)); meany=mean(points(:,2)); meanz=mean(points(:,3));
%             %Ddis1=sqrt((meanx-cellcenter(1,1))^2+(meany-cellcenter(1,2))^2);
%             msd2 = MSD_value(frameNum,points(:,1:2), MSD_length);%*0.2667*0.2667;
%             [DF2,AlphaF2,RsquareF2,DL2,AlphaL2,RsquareL2]=MSDtrack2dF4L(msd2,4,33,TimeResolution);
%             %[D2,Alpha2,Rsquare2]=MSDtrack2d(msd2,3,alphasize,0.03);
%             msd3 = MSD_value3d(frameNum,points(:,1:3), MSD_length);%*0.2667*0.2667*0.18;
%             [DF3,AlphaF3,RsquareF3,DL3,AlphaL3,RsquareL3]=MSDtrack3dF4L(msd3,4,33,TimeResolution);
%             %[D3,Alpha3,Rsquare3]=MSDtrack3d(msd3,3,alphasize,0.03);       
%             msd1 = MSD_value1d(frameNum,points(:,3), MSD_length);
%             [DF1,AlphaF1,RsquareF1,DL1,AlphaL1,RsquareL1]=MSDtrack1dF4L(msd1,4,33,TimeResolution);
%             %[D1,Alpha1,Rsquare1]=MSDtrack1d(msd1,3,alphasize,0.03);
%             Da1=[meanx,meany,meanz,DF3,AlphaF3,DL3,AlphaL3,DF2,AlphaF2,DL2,AlphaL2,DF1,AlphaF1,DL1,AlphaL1];
%             %Ddis=[Ddis;Ddis1];
%             DaF4L=[DaF4L;Da1];
%             end
%     end
    Da2=[];
    DaxyzF4L=[];                                                                                                                                                                                                                                                                                                                                                                                                                           
    for itrack = 1:length(Track)
            frameNum=Track(itrack,1).frameNum;
            points=Track(itrack,1).points;
            %intensities=Track(itrack,1).intensities;
            if length(points(:,1))>limit_length
            meanx=mean(points(:,1)); meany=mean(points(:,2)); meanz=mean(points(:,3));
            msdx = MSD_value1d(frameNum,points(:,1), MSD_length);%*0.2667*0.2667;
            %[Dx,Alphax,Rsquarex]=MSDtrack1d(msdx,3,alphasize,timelag);
            [DxF,AlphaxF,DxL,AlphaxL]=MSDtrack1dF4LDAjc(msdx,7,33,TimeResolution);
            msdy = MSD_value1d(frameNum,points(:,2), MSD_length);%*0.2667*0.2667*0.18;
            [DyF,AlphayF,DyL,AlphayL]=MSDtrack1dF4LDAjc(msdy,7,33,TimeResolution);
            %[Dy,Alphay,Rsquarey]=MSDtrack1d(msdy,3,alphasize,timelag);       
            msdz = MSD_value1d(frameNum,points(:,3), MSD_length);
            [DzF,AlphazF,DzL,AlphazL]=MSDtrack1dF4LDAjc(msdz,7,33,TimeResolution);
            %if DxF<0.5&DxF>0&DxL<0.5&DxL>0&DyF<0.5&DyF>0&DyL<0.5&DyL>0&DzF<0.5&DzF>0&DzL<0.5&DzL>0&...
                   % AlphaxF<2&AlphaxF>0&AlphaxL<2&AlphaxL>0&AlphayF<2&AlphayF>0&AlphayL<2&AlphayL>0&AlphazF<2&AlphazF>0&AlphazL<2&AlphazL>0
            %if   AlphaxF>0.1&AlphayF>0.1&AlphazF>0.1
                    
            %[Dz,Alphaz,Rsquarez]=MSDtrack1d(msdz,3,alphasize,timelag);
            Da2=[DxF,AlphaxF,DxL,AlphaxL,DyF,AlphayF,DyL,AlphayL,DzF,AlphazF,DzL,AlphazL];
            DaxyzF4L=[DaxyzF4L;Da2];
            %end           
            end
    end
    else
        DaxyzF4L=[];DaF4L=[];
    end
    %ViscoCellData.DaF4L=DaF4L;
    ViscoCellData.DaxyzF4L=DaxyzF4L;
    ViscoCellData1(iSample,1).ViscoCellData=ViscoCellData;
    clear ViscoCellData;
end
%%disp(strcat('down-Next!'));
%%scatter(Da(:,7),Da(:,10));