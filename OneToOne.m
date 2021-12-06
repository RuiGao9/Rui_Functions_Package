function [Rsquare_org,RMSE] = OneToOne(x,y,size_scatter,...
    color_scatter,color_line,XLABEL,YLABEL)
%ONETOONE Summary of this function goes here
%   Detailed explanation goes here
%   Parameters:
%   x and y: two vectors
%   size_scatter: the size of the scatters shown on the figure
%   color_scatter: scatter color
%   color_line: the color of the trend line (linear regression)
%   XLABEL and YLABEL: the label shown on the figure
%   Outputs:
%   Rsquare_org: R-square value
%   RMSE: Root mean square error

mdl = fitlm(x,y); 
Rsquare_org = mdl.Rsquared.Ordinary;
RMSE = mdl.RMSE;

b1 = mdl.Coefficients.Estimate(2);
inter = mdl.Coefficients.Estimate(1);
yCal = b1*x+inter;

figure
scatterhist(x,y,'Color','b','Kernel','on','Marker','+'); hold on
plot(x,yCal,color_line); hold off
xlabel(XLABEL); ylabel(YLABEL);

text(max(x)*0.95,max(y)*0.95,['RMSE:',num2str(round(RMSE,2)),"R-square:",num2str(round(Rsquare_org,2))],'FontSize',12)

print('-dpng');
end

