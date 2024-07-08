function fig=selectionHelperPlotBETA(beta)
reps=size(beta,4);
dynLvl=size(beta,1);
totalPlots=dynLvl*reps;
span=1:reps:totalPlots;
fig=figure('WindowState','maximized');
T=tiledlayout(fig,dynLvl,reps);
T.XLabel.String="Components";
T.YLabel.String="Beta";
T.XLabel.FontSize=20;
T.YLabel.FontSize=20;
oldFolder=cd([toolboxdir('stats'),'\stats']);
warning('off')
for i=1:dynLvl
    for j=1:reps
        temp=permute(beta(i,:,:,j),[3,2,1]);
        comps=size(temp,2);


        t1=nexttile(T,span(i)+j-1);
        boxplot(temp.^2)
        %         boxplot(temp)
        hold on
        yline(0)
        text(1.2:1:comps+0.2,0.90*ones(1,comps),cellstr(num2str(var(temp)')))
        t1.XTickLabel=[];
        t1.YLim=[0 1];
        %         t1.YLim=[-1 1];
        t1.YLabel.String=num2str(i-1);
        t1.YLabel.FontSize=18;
        box on
    end
end
warning('on')
cd(oldFolder)