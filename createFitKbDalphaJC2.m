function [fitresult, gof] =createFitKbDalphaJC2(x, y)
%CREATEFIT(TF,MSDTTXF)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : tF
%      Y Output: MSDttxF
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 27-Dec-2021 16:55:38 自动生成


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
ft = fittype( '2*d*(x^a)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Algorithm = 'Levenberg-Marquardt';
opts.Display = 'Off';
opts.StartPoint = [0.2 0.01];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% % Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'MSDttxF vs. tF', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% xlabel tF
% ylabel MSDttxF
% grid on


