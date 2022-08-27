function [fitresult, gof] = createFitradius(Y, X, I,r)
% Y，X分别是meshgrid给出的X，Y，对应关系相反
%CREATEFIT(Y,X,I)

%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : Y
%      Y Input : X
%      Z Output: I
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 02-Feb-2018 22:28:27 自动生成


%% Fit: 'untitled fit 1'.
[xData, yData, zData] = prepareSurfaceData( Y, X, I );

% Set up fittype and options.
ft = fittype( 'a + b*exp(-d*((x-m).^2+(y-n).^2))+c*exp(-e*((((x-m).^2+(y-n).^2)^0.5-p).^2))', 'independent', {'x', 'y'}, 'dependent', 'z' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.DiffMaxChange = 0.2;
opts.DiffMinChange = 1e-09;
opts.Display = 'Off';
opts.Lower = [-5 -5 -5 -5 -5 -7 -7 0.5];
opts.Robust = 'Bisquare';
opts.StartPoint = [0.1 0.1 0.1 0.1 0.1 0.1 0.1 r];
opts.Upper = [5 5 5 5 5 7 7 24.1];

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft, opts );




