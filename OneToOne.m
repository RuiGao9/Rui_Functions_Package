function [Rsquare_org,RMSE] = OneToOne(x,y,...
    color_line,XLABEL,YLABEL,...
    modelspec)
%ONETOONE Summary of this function goes here
%   Detailed explanation goes here
%   Parameters:
%   x and y: two vectors
%   color_line: the color of the trend line (linear regression)
%   XLABEL and YLABEL: the label shown on the figure
%   modelspec: 2 types are available now.
%               "linear": default in this function. Model contains an intercept and linear term for each predictor.
%               "purequadratic": Model contains an intercept term and linear and squared terms for each predictor.
%   Outputs:
%   Rsquare_org: R-square value
%   RMSE: Root mean square error

if modelspec == "linear"
    mdl = fitlm(x,y,"linear");
    
    Rsquare_org = mdl.Rsquared.Ordinary;
    RMSE = mdl.RMSE;
    b1 = mdl.Coefficients.Estimate(2);
    inter = mdl.Coefficients.Estimate(1);
    yCal = b1.*x+inter;
elseif modelspec == "purequadratic"
    mdl = fitlm(x,y,"purequadratic")
    
    Rsquare_org = mdl.Rsquared.Ordinary;
    RMSE = mdl.RMSE;
    b2 = mdl.Coefficients.Estimate(3);
    b1 = mdl.Coefficients.Estimate(2);
    inter = mdl.Coefficients.Estimate(1);
    yCal = b2.*x.^2+b1.*x+inter;
end

table_plot = [x,yCal];
table_plot = sortrows(table_plot,1);
figure
scatterhist(x,y,'Color','b','Kernel','on','Marker','+'); hold on
plot(table_plot(:,1),table_plot(:,2),color_line); hold off
xlabel(XLABEL); ylabel(YLABEL);

% text(max(x)*0.95,max(y)*0.95,['RMSE:',num2str(round(RMSE,2)),"R-square:",num2str(round(Rsquare_org,2))],'FontSize',12)
yL = ylim;
text(max(x),0.95*yL(2),...
    ['RMSE:',num2str(round(RMSE,2)),"R-square:",num2str(round(Rsquare_org,2))],'FontSize',12)

print('-dpng');
end

