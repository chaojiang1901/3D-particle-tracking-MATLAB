function [DF,AlphaF,DL,AlphaL]=MSDtrack1dfit(msd,D_MSDSize,Alpha_MSDSize,TimeResolution);
    %D_MSDSize=3;     
    dt1=(1:100)*TimeResolution; dt=dt1';
    %ft_ = fittype('poly1');
    x=dt(1:D_MSDSize);y=msd(1:D_MSDSize);
    [fitresultF, gofF] = createFitKbDalpha(x, y); %JC1 MSD=d*t^a+c;
    %[cf_,gof] = fit(dt(1:D_MSDSize),msd(1:D_MSDSize),ft_);  %多项式拟合
    DF = fitresultF.d; %*PixelSize*PixelSize;
    AlphaF=fitresultF.a;
    %Alpha_MSDSize=10;
%     lg_msd = log(msd(1:D_MSDSize));
%     lg_t=log(dt(1:D_MSDSize));
%     ok_ = isfinite(lg_t)&isfinite(lg_msd);
%     ft_ = fittype('poly1');
%     [cf_,gof] = fit(lg_t(ok_),lg_msd(ok_),ft_);
%     AlphaF = cf_.p1;
%     RsquareF = gof.rsquare;
     x1=dt(D_MSDSize+1:Alpha_MSDSize);
     y1=msd(D_MSDSize+1:Alpha_MSDSize);
    [fitresultL, gofL] = createFitKbDalpha(x1, y1);
    DL = fitresultL.d; %*PixelSize*PixelSize;
    AlphaL=fitresultL.a;
    %dt1=(1:100)*TimeResolution; dt=dt1';
%     ft_ = fittype('poly1');
%     [cf_,gof] = fit(dt(D_MSDSize+1:Alpha_MSDSize),msd(D_MSDSize+1:Alpha_MSDSize),ft_);  %多项式拟合
%     DL = cf_.p1/2; %*PixelSize*PixelSize;
%     %Alpha_MSDSize=10;
%     lg_msd = log(msd(D_MSDSize+1:Alpha_MSDSize));
%     lg_t=log(dt(D_MSDSize+1:Alpha_MSDSize));
%     ok_ = isfinite(lg_t)&isfinite(lg_msd);
%     ft_ = fittype('poly1');
%     [cf_,gof] = fit(lg_t(ok_),lg_msd(ok_),ft_);
%     AlphaL= cf_.p1;
%     RsquareL = gof.rsquare;