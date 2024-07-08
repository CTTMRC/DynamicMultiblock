load('.\Data\TabularResults.mat',"fullResultsR2")
rowNames=fullResultsR2.Properties.RowNames;
rowNames=regexp(rowNames,"_",'split');
rowNames=cellfun(@(x)x{:,1},rowNames,'UniformOutput',false);
rowNames=split(rowNames,"-");
varNames=fullResultsR2.Properties.VariableNames';
groupNames=[repmat(rowNames,7,1),repelem(varNames,24,1)];
dynamicStrenght=double(categorical((groupNames(:,1))));
Contribution=double(categorical((groupNames(:,2))));
spectraComplexity=double(categorical((groupNames(:,3))));
modelName=double(categorical((groupNames(:,4))));
dynamicStrenght=categorical(dynamicStrenght,unique(dynamicStrenght,'stable'),...
     ["[0;0]","[1;1]","[0;2]","[2;0]"]);
Contribution=categorical(Contribution,unique(Contribution,'stable'),...
     ["[0.5;0.5]","[0.8;0.2]","[0.2;0.8]"]);
spectraComplexity=categorical(spectraComplexity,unique(spectraComplexity,'stable'),...
 ["[1]","[2345]"]);
modelName=categorical(modelName,unique(modelName,'stable'),...
 varNames);
levels=reshape(fullResultsR2.Variables,[],1);
%%
[~,T,STATS] = anovan(levels,{dynamicStrenght,Contribution,spectraComplexity,modelName},"model",3,...
'varnames',{'DynOrder','blockContribution','spectraComplexity', 'ModelType' });
rowNames=T(2:end,1);
tvarNames=T(1,2:end);
T=T(2:end,2:end);
T{end,4}=0;
T{end-1,5}=0;
T{end,5}=0;
T{end-1,6}=0;
T{end,6}=0;
T=cell2table(T);
T.Properties.VariableNames=tvarNames;
T.Properties.RowNames=rowNames;
T.('Sum Sq.')=round(T.('Sum Sq.'),3);
T.('Mean Sq.')=round(T.('Mean Sq.'),5);
T.('F')=round(T.('F'),3);
T.('Prob>F')=round(T.('Prob>F'),3);
T.('Explained Variance')=round((T.('Sum Sq.')./T.('Sum Sq.')(end))*100,2);
T.('Singular?')=[];
writetable(T,"./Data/anovaResults.xlsx",'WriteRowNames',true,'WriteVariableNames',true);
figure();
[~,~,f]=multcompare(STATS,'Dimension',4,'CType','lsd');
plotData=readMultcomparePlot(f);
[~,I]=sort(plotData.m,'ascend');
performanceOrder=double(categorical(plotData.ModelType(I)));

[~,~,f]=multcompare(STATS,'Dimension',[1,2,4],'CType','lsd');
plotData=readMultcomparePlot(f);
figure();
S=cellstr(["v","s","o","d","^","x"]);
C=cellstr(["r","b","#017d03","#d47902"]);
MSz=8;
f1=figure();
T1=tiledlayout(f1,3,1);
T1.TileSpacing = 'tight';
T1.Padding = 'compact';
t1=nexttile(1);
tmp = (plotData(plotData.blockContribution=="[0.5;0.5]",:));%,:);%
catMethods=categorical(tmp.ModelType,varNames);
offset=[-0.15,-0.05,0.05,0.15];%[-0.2,-0.1,0.1,0.2];
UL=nan(4,7);
h=gscatter(double(catMethods),tmp.ll,tmp.DynOrder); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = '_';
    h(i).MarkerSize= MSz*2;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = C{i} ;
    h(i).Color = 'k';%C(i);
    UL(i,:)=h(i).YData;
end
hold on
h=gscatter(double(catMethods),tmp.ul,tmp.DynOrder); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = '_';
    h(i).MarkerSize= MSz*2;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = C{i} ;
    h(i).Color = 'k';%C(i);
    for j=1:size(h(1).XData,2)
        line([h(i).XData(j),h(i).XData(j)],[UL(i,j),h(i).YData(j)],'Color','k')
    end
end
h=gscatter(double(catMethods),tmp.m,tmp.DynOrder); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = S{i};
    h(i).MarkerSize= MSz;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = 'k';%C(i) ;
    h(i).Color = 'k';%C(i);
end
t1.Title.String='blockContribution=="[0.5;0.5]"';
t1.Title.FontSize=26;
t1.YLim=[0 1];
t1.XLim=[0.5 7.5];
t1.XTick=1:7;
t1.XTickLabel=categories(catMethods);
t1.XTickLabelRotation=30;
t1.FontSize=18;
box on;

t3=nexttile(2);
tmp = (plotData(plotData.blockContribution=="[0.8;0.2]",:));%,:);%
catMethods=categorical(tmp.ModelType,varNames);
offset=[-0.15,-0.05,0.05,0.15];%[-0.2,-0.1,0.1,0.2];
UL=nan(4,7);
h=gscatter(double(catMethods),tmp.ll,tmp.DynOrder); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = '_';
    h(i).MarkerSize= MSz*2;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = C{i} ;
    h(i).Color = 'k';%C(i);
    UL(i,:)=h(i).YData;
end
hold on
h=gscatter(double(catMethods),tmp.ul,tmp.DynOrder); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = '_';
    h(i).MarkerSize= MSz*2;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = C{i} ;
    h(i).Color = 'k';%C(i);
    for j=1:size(h(1).XData,2)
        line([h(i).XData(j),h(i).XData(j)],[UL(i,j),h(i).YData(j)],'Color','k')
    end
end
h=gscatter(double(catMethods),tmp.m,tmp.DynOrder); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = S{i};
    h(i).MarkerSize= MSz;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = 'k';%C(i) ;
    h(i).Color = 'k';%C(i);
end
t3.Title.String='blockContribution=="[0.8;0.2]"';
t3.Title.FontSize=26;
t3.YLim=[0 1];
t3.XLim=[0.5 7.5];
t3.XTick=1:7;
t3.XTickLabel=categories(catMethods);
t3.XTickLabelRotation=30;
t3.FontSize=18;
box on;

t5=nexttile(3);
tmp = (plotData(plotData.blockContribution=="[0.2;0.8]",:));%,:);%
catMethods=categorical(tmp.ModelType,varNames);
offset=[-0.15,-0.05,0.05,0.15];%[-0.2,-0.1,0.1,0.2];
UL=nan(4,7);
h=gscatter(double(catMethods),tmp.ll,tmp.DynOrder); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = '_';
    h(i).MarkerSize= MSz*2;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = C{i} ;
    h(i).Color = 'k';%C(i);
    UL(i,:)=h(i).YData;
end
hold on
h=gscatter(double(catMethods),tmp.ul,tmp.DynOrder); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = '_';
    h(i).MarkerSize= MSz*2;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = C{i} ;
    h(i).Color = 'k';%C(i);
    for j=1:size(h(1).XData,2)
        line([h(i).XData(j),h(i).XData(j)],[UL(i,j),h(i).YData(j)],'Color','k')
    end
end
h=gscatter(double(catMethods),tmp.m,tmp.DynOrder); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = S{i};
    h(i).MarkerSize= MSz;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = 'k';%C(i) ;
    h(i).Color = 'k';%C(i);
end
t5.Title.String='blockContribution=="[0.2;0.8]"';
t5.Title.FontSize=26;
t5.YLim=[0 1];
t5.XLim=[0.5 7.5];
t5.XTick=1:7;
t5.XTickLabel=categories(catMethods);
t5.XTickLabelRotation=30;
t5.FontSize=18;
box on;
leg=legend(h,'Orientation', 'vertical');
leg.Layout.Tile = 'east';
leg.Interpreter = 'tex';
leg.Title.String="DynOrder";
box on;
T1.XLabel.String="Method";T1.XLabel.FontWeight="bold";T1.XLabel.FontSize=26;
T1.YLabel.String="R^2_P"+newline;T1.YLabel.FontWeight="bold";T1.YLabel.FontSize=26;
%%
[P,T,STATS] = anovan(levels,{dynamicStrenght,Contribution,spectraComplexity,modelName},"model",3,...
'varnames',{'DynOrder','blockContribution','spectraComplexity', 'ModelType' });
figure();
[~,~,f]=multcompare(STATS,'Dimension',4,'CType','lsd');
plotData=readMultcomparePlot(f);
[~,I]=sortrows(plotData,"m",'ascend');
performanceOrder=double(categorical(plotData.ModelType(I)));
uniqueMethods=unique(plotData.ModelType(I),'stable');

[~,~,f]=multcompare(STATS,'Dimension',[1,2,4],'CType','lsd');
plotData=readMultcomparePlot(f);
[~,I]=sortrows(plotData,"m",'ascend');
figure();
S=cellstr(["v","s","o","d","^","x"]);
C=cellstr(["r","b","#017d03","#d47902"]);
MSz=8;
f1=figure();
T1=tiledlayout(f1,4,1);
T1.TileSpacing = 'tight';
T1.Padding = 'compact';
t1=nexttile(1);
tmp = (plotData(plotData.DynOrder=="[0;0]",:));%,:);%
catMethods=categorical(tmp.ModelType,varNames);
offset=[-0.15,-0.05,0.05,0.15];%[-0.1,0,0.1];
UL=nan(3,7);
h=gscatter(double(catMethods),tmp.ll,tmp.blockContribution); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = '_';
    h(i).MarkerSize= MSz*2;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = C{i} ;
    h(i).Color = 'k';%C(i);
    UL(i,:)=h(i).YData;
end
hold on
h=gscatter(double(catMethods),tmp.ul,tmp.blockContribution); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = '_';
    h(i).MarkerSize= MSz*2;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = C{i} ;
    h(i).Color = 'k';%C(i);
    for j=1:size(h(1).XData,2)
        line([h(i).XData(j),h(i).XData(j)],[UL(i,j),h(i).YData(j)],'Color','k')
    end
end
h=gscatter(double(catMethods),tmp.m,tmp.blockContribution); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = S{i};
    h(i).MarkerSize= MSz;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = 'k';%C(i) ;
    h(i).Color = 'k';%C(i);
end
t1.Title.String='DynOrder=="[0;0]"';
t1.Title.FontSize=26;
t1.YLim=[0 1];
t1.XLim=[0.5 7.5];
t1.XTick=1:7;
t1.XTickLabel=categories(catMethods);
t1.XTickLabelRotation=30;
t1.FontSize=18;
box on;

t3=nexttile(2);
tmp = (plotData(plotData.DynOrder=="[1;1]",:));%,:);%
catMethods=categorical(tmp.ModelType,varNames);
offset=[-0.15,-0.05,0.05,0.15];%[-0.1,0,0.1];
UL=nan(3,7);
h=gscatter(double(catMethods),tmp.ll,tmp.blockContribution); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = '_';
    h(i).MarkerSize= MSz*2;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = C{i} ;
    h(i).Color = 'k';%C(i);
    UL(i,:)=h(i).YData;
end
hold on
h=gscatter(double(catMethods),tmp.ul,tmp.blockContribution); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = '_';
    h(i).MarkerSize= MSz*2;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = C{i} ;
    h(i).Color = 'k';%C(i);
    for j=1:size(h(1).XData,2)
        line([h(i).XData(j),h(i).XData(j)],[UL(i,j),h(i).YData(j)],'Color','k')
    end
end
h=gscatter(double(catMethods),tmp.m,tmp.blockContribution); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = S{i};
    h(i).MarkerSize= MSz;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = 'k';%C(i) ;
    h(i).Color = 'k';%C(i);
end
t3.Title.String='DynOrder=="[1;1]"';
t3.Title.FontSize=26;
t3.YLim=[0 1];
t3.XLim=[0.5 7.5];
t3.XTick=1:7;
t3.XTickLabel=categories(catMethods);
t3.XTickLabelRotation=30;
t3.FontSize=18;
box on;

t5=nexttile(3);
tmp = (plotData(plotData.DynOrder=="[2;0]",:));%,:);%
catMethods=categorical(tmp.ModelType,varNames);
offset=[-0.15,-0.05,0.05,0.15];%[-0.1,0,0.1];
UL=nan(3,7);
h=gscatter(double(catMethods),tmp.ll,tmp.blockContribution); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = '_';
    h(i).MarkerSize= MSz*2;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = C{i} ;
    h(i).Color = 'k';%C(i);
    UL(i,:)=h(i).YData;
end
hold on
h=gscatter(double(catMethods),tmp.ul,tmp.blockContribution); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = '_';
    h(i).MarkerSize= MSz*2;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = C{i} ;
    h(i).Color = 'k';%C(i);
    for j=1:size(h(1).XData,2)
        line([h(i).XData(j),h(i).XData(j)],[UL(i,j),h(i).YData(j)],'Color','k')
    end
end
h=gscatter(double(catMethods),tmp.m,tmp.blockContribution); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = S{i};
    h(i).MarkerSize= MSz;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = 'k';%C(i) ;
    h(i).Color = 'k';%C(i);
end
t5.Title.String='DynOrder=="[2;0]"';
t5.Title.FontSize=26;
t5.YLim=[0 1];
t5.XLim=[0.5 7.5];
t5.XTick=1:7;
t5.XTickLabel=categories(catMethods);
t5.XTickLabelRotation=30;
t5.FontSize=18;
box on;

t7=nexttile(4);
tmp = (plotData(plotData.DynOrder=="[0;2]",:));%,:);%
catMethods=categorical(tmp.ModelType,varNames);
offset=[-0.15,-0.05,0.05,0.15];%[-0.1,0,0.1];
UL=nan(3,7);
h=gscatter(double(catMethods),tmp.ll,tmp.blockContribution); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = '_';
    h(i).MarkerSize= MSz*2;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = C{i} ;
    h(i).Color = 'k';%C(i);
    UL(i,:)=h(i).YData;
end
hold on
h=gscatter(double(catMethods),tmp.ul,tmp.blockContribution); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = '_';
    h(i).MarkerSize= MSz*2;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = C{i} ;
    h(i).Color = 'k';%C(i);
    for j=1:size(h(1).XData,2)
        line([h(i).XData(j),h(i).XData(j)],[UL(i,j),h(i).YData(j)],'Color','k')
    end
end
h=gscatter(double(catMethods),tmp.m,tmp.blockContribution); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = S{i};
    h(i).MarkerSize= MSz;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = 'k';%C(i) ;
    h(i).Color = 'k';%C(i);
end
t7.Title.String='DynOrder=="[0;2]"';
t7.Title.FontSize=26;
t7.YLim=[0 1];
t7.XLim=[0.5 7.5];
t7.XTick=1:7;
t7.XTickLabel=categories(catMethods);
t7.XTickLabelRotation=30;
t7.FontSize=18;
box on;

leg=legend(h,'Orientation', 'vertical');
leg.Layout.Tile = 'east';
leg.Interpreter = 'tex';
leg.Title.String="blockContribution";
box on;
T1.XLabel.String="Method";T1.XLabel.FontWeight="bold";T1.XLabel.FontSize=26;
T1.YLabel.String="R^2_P"+newline;T1.YLabel.FontWeight="bold";T1.YLabel.FontSize=26;



%%
load('.\Data\TabularResults.mat',"fullResultsR2")
fullResultsR2.PLS=[];
fullResultsR2.SOPLS=[];
rowNames=fullResultsR2.Properties.RowNames;
rowNames=regexp(rowNames,"_",'split');
rowNames=cellfun(@(x)x{:,1},rowNames,'UniformOutput',false);
rowNames=split(rowNames,"-");
varNames=fullResultsR2.Properties.VariableNames';
groupNames=[repmat(rowNames,5,1),repelem(varNames,24,1)];
dynamicStrenght=double(categorical((groupNames(:,1))));
Contribution=double(categorical((groupNames(:,2))));
spectraComplexity=double(categorical((groupNames(:,3))));
modelName=double(categorical((groupNames(:,4))));
dynamicStrenght=categorical(dynamicStrenght,unique(dynamicStrenght,'stable'),...
     ["[0;0]","[1;1]","[0;2]","[2;0]"]);
Contribution=categorical(Contribution,unique(Contribution,'stable'),...
     ["[0.5;0.5]","[0.8;0.2]","[0.2;0.8]"]);
spectraComplexity=categorical(spectraComplexity,unique(spectraComplexity,'stable'),...
 ["[1]","[2345]"]);
modelName=categorical(modelName,unique(modelName,'stable'),...
 varNames);
levels=reshape(fullResultsR2.Variables,[],1);
%%
[~,T,STATS] = anovan(levels,{dynamicStrenght,Contribution,spectraComplexity,modelName},"model",3,...
'varnames',{'DynOrder','blockContribution','spectraComplexity', 'ModelType' });
rowNames=T(2:end,1);
tvarNames=T(1,2:end);
T=T(2:end,2:end);
T{end,4}=0;
T{end-1,5}=0;
T{end,5}=0;
T{end-1,6}=0;
T{end,6}=0;
T=cell2table(T);
T.Properties.VariableNames=tvarNames;
T.Properties.RowNames=rowNames;
T.('Sum Sq.')=round(T.('Sum Sq.'),3);
T.('Mean Sq.')=round(T.('Mean Sq.'),5);
T.('F')=round(T.('F'),3);
T.('Prob>F')=round(T.('Prob>F'),3);
T.('Explained Variance')=round((T.('Sum Sq.')./T.('Sum Sq.')(end))*100,2);
T.('Singular?')=[];
writetable(T,"./Data/anovaResultsDO.xlsx",'WriteRowNames',true,'WriteVariableNames',true);
figure();
[~,~,f]=multcompare(STATS,'Dimension',4,'CType','lsd');
plotData=readMultcomparePlot(f);
[~,I]=sortrows(plotData,"m",'ascend');
performanceOrder=double(categorical(plotData.ModelType(I)));

[~,~,f]=multcompare(STATS,'Dimension',[1,2,4],'CType','lsd');
plotData=readMultcomparePlot(f);
figure();
S=cellstr(["v","s","o","d","^","x"]);
C=cellstr(["r","b","#017d03","#d47902"]);
MSz=8;
f1=figure();
T1=tiledlayout(f1,3,1);
T1.TileSpacing = 'tight';
T1.Padding = 'compact';
t1=nexttile(1);
tmp = (plotData(plotData.blockContribution=="[0.5;0.5]",:));%,:);%
catMethods=categorical(tmp.ModelType,varNames);
offset=[-0.15,-0.05,0.05,0.15];%[-0.2,-0.1,0.1,0.2];
UL=nan(4,5);
h=gscatter(double(catMethods),tmp.ll,tmp.DynOrder); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = '_';
    h(i).MarkerSize= MSz*2;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = C{i} ;
    h(i).Color = 'k';%C(i);
    UL(i,:)=h(i).YData;
end
hold on
h=gscatter(double(catMethods),tmp.ul,tmp.DynOrder); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = '_';
    h(i).MarkerSize= MSz*2;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = C{i} ;
    h(i).Color = 'k';%C(i);
    for j=1:size(h(1).XData,2)
        line([h(i).XData(j),h(i).XData(j)],[UL(i,j),h(i).YData(j)],'Color','k')
    end
end
h=gscatter(double(catMethods),tmp.m,tmp.DynOrder); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = S{i};
    h(i).MarkerSize= MSz;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = 'k';%C(i) ;
    h(i).Color = 'k';%C(i);
end
t1.Title.String='blockContribution=="[0.5;0.5]"';
t1.Title.FontSize=26;
t1.YLim=[0 1];
t1.XLim=[0.5 5.5];
t1.XTick=1:5;
t1.XTickLabel=categories(catMethods);
t1.XTickLabelRotation=30;
t1.FontSize=18;
box on;

t3=nexttile(2);
tmp = (plotData(plotData.blockContribution=="[0.8;0.2]",:));%,:);%
catMethods=categorical(tmp.ModelType,varNames);
offset=[-0.15,-0.05,0.05,0.15];%[-0.2,-0.1,0.1,0.2];
UL=nan(4,5);
h=gscatter(double(catMethods),tmp.ll,tmp.DynOrder); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = '_';
    h(i).MarkerSize= MSz*2;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = C{i} ;
    h(i).Color = 'k';%C(i);
    UL(i,:)=h(i).YData;
end
hold on
h=gscatter(double(catMethods),tmp.ul,tmp.DynOrder); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = '_';
    h(i).MarkerSize= MSz*2;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = C{i} ;
    h(i).Color = 'k';%C(i);
    for j=1:size(h(1).XData,2)
        line([h(i).XData(j),h(i).XData(j)],[UL(i,j),h(i).YData(j)],'Color','k')
    end
end
h=gscatter(double(catMethods),tmp.m,tmp.DynOrder); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = S{i};
    h(i).MarkerSize= MSz;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = 'k';%C(i) ;
    h(i).Color = 'k';%C(i);
end
t3.Title.String='blockContribution=="[0.8;0.2]"';
t3.Title.FontSize=26;
t3.YLim=[0 1];
t3.XLim=[0.5 5.5];
t3.XTick=1:5;
t3.XTickLabel=categories(catMethods);
t3.XTickLabelRotation=30;
t3.FontSize=18;
box on;

t5=nexttile(3);
tmp = (plotData(plotData.blockContribution=="[0.2;0.8]",:));%,:);%
catMethods=categorical(tmp.ModelType,varNames);
offset=[-0.15,-0.05,0.05,0.15];%[-0.2,-0.1,0.1,0.2];
UL=nan(4,5);
h=gscatter(double(catMethods),tmp.ll,tmp.DynOrder); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = '_';
    h(i).MarkerSize= MSz*2;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = C{i} ;
    h(i).Color = 'k';%C(i);
    UL(i,:)=h(i).YData;
end
hold on
h=gscatter(double(catMethods),tmp.ul,tmp.DynOrder); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = '_';
    h(i).MarkerSize= MSz*2;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = C{i} ;
    h(i).Color = 'k';%C(i);
    for j=1:size(h(1).XData,2)
        line([h(i).XData(j),h(i).XData(j)],[UL(i,j),h(i).YData(j)],'Color','k')
    end
end
h=gscatter(double(catMethods),tmp.m,tmp.DynOrder); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = S{i};
    h(i).MarkerSize= MSz;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = 'k';%C(i) ;
    h(i).Color = 'k';%C(i);
end
t5.Title.String='blockContribution=="[0.2;0.8]"';
t5.Title.FontSize=26;
t5.YLim=[0 1];
t5.XLim=[0.5 5.5];
t5.XTick=1:5;
t5.XTickLabel=categories(catMethods);
t5.XTickLabelRotation=30;
t5.FontSize=18;
box on;
leg=legend(h,'Orientation', 'vertical');
leg.Layout.Tile = 'east';
leg.Interpreter = 'tex';
leg.Title.String="DynOrder";
box on;
T1.XLabel.String="Method";T1.XLabel.FontWeight="bold";T1.XLabel.FontSize=26;
T1.YLabel.String="R^2_P"+newline;T1.YLabel.FontWeight="bold";T1.YLabel.FontSize=26;
%%
[P,T,STATS] = anovan(levels,{dynamicStrenght,Contribution,spectraComplexity,modelName},"model",3,...
'varnames',{'DynOrder','blockContribution','spectraComplexity', 'ModelType' });
figure();
[~,~,f]=multcompare(STATS,'Dimension',4,'CType','lsd');
plotData=readMultcomparePlot(f);
[~,I]=sortrows(plotData,"m",'ascend');
performanceOrder=double(categorical(plotData.ModelType(I)));
uniqueMethods=unique(plotData.ModelType(I),'stable');

[~,~,f]=multcompare(STATS,'Dimension',[1,2,4],'CType','lsd');
plotData=readMultcomparePlot(f);
[~,I]=sortrows(plotData,"m",'ascend');
figure();
S=cellstr(["v","s","o","d","^","x"]);
C=cellstr(["r","b","#017d03","#d47902"]);
MSz=8;
f1=figure();
T1=tiledlayout(f1,4,1);
T1.TileSpacing = 'tight';
T1.Padding = 'compact';
t1=nexttile(1);
tmp = (plotData(plotData.DynOrder=="[0;0]",:));%,:);%
catMethods=categorical(tmp.ModelType,varNames);
offset=[-0.15,-0.05,0.05,0.15];%[-0.1,0,0.1];
UL=nan(3,5);
h=gscatter(double(catMethods),tmp.ll,tmp.blockContribution); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = '_';
    h(i).MarkerSize= MSz*2;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = C{i} ;
    h(i).Color = 'k';%C(i);
    UL(i,:)=h(i).YData;
end
hold on
h=gscatter(double(catMethods),tmp.ul,tmp.blockContribution); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = '_';
    h(i).MarkerSize= MSz*2;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = C{i} ;
    h(i).Color = 'k';%C(i);
    for j=1:size(h(1).XData,2)
        line([h(i).XData(j),h(i).XData(j)],[UL(i,j),h(i).YData(j)],'Color','k')
    end
end
h=gscatter(double(catMethods),tmp.m,tmp.blockContribution); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = S{i};
    h(i).MarkerSize= MSz;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = 'k';%C(i) ;
    h(i).Color = 'k';%C(i);
end
t1.Title.String='DynOrder=="[0;0]"';
t1.Title.FontSize=26;
t1.YLim=[0 1];
t1.XLim=[0.5 5.5];
t1.XTick=1:5;
t1.XTickLabel=categories(catMethods);
t1.XTickLabelRotation=30;
t1.FontSize=18;
box on;

t3=nexttile(2);
tmp = (plotData(plotData.DynOrder=="[1;1]",:));%,:);%
catMethods=categorical(tmp.ModelType,varNames);
offset=[-0.15,-0.05,0.05,0.15];%[-0.1,0,0.1];
UL=nan(3,5);
h=gscatter(double(catMethods),tmp.ll,tmp.blockContribution); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = '_';
    h(i).MarkerSize= MSz*2;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = C{i} ;
    h(i).Color = 'k';%C(i);
    UL(i,:)=h(i).YData;
end
hold on
h=gscatter(double(catMethods),tmp.ul,tmp.blockContribution); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = '_';
    h(i).MarkerSize= MSz*2;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = C{i} ;
    h(i).Color = 'k';%C(i);
    for j=1:size(h(1).XData,2)
        line([h(i).XData(j),h(i).XData(j)],[UL(i,j),h(i).YData(j)],'Color','k')
    end
end
h=gscatter(double(catMethods),tmp.m,tmp.blockContribution); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = S{i};
    h(i).MarkerSize= MSz;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = 'k';%C(i) ;
    h(i).Color = 'k';%C(i);
end
t3.Title.String='DynOrder=="[1;1]"';
t3.Title.FontSize=26;
t3.YLim=[0 1];
t3.XLim=[0.5 5.5];
t3.XTick=1:5;
t3.XTickLabel=categories(catMethods);
t3.XTickLabelRotation=30;
t3.FontSize=18;
box on;

t5=nexttile(3);
tmp = (plotData(plotData.DynOrder=="[2;0]",:));%,:);%
catMethods=categorical(tmp.ModelType,varNames);
offset=[-0.15,-0.05,0.05,0.15];%[-0.1,0,0.1];
UL=nan(3,5);
h=gscatter(double(catMethods),tmp.ll,tmp.blockContribution); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = '_';
    h(i).MarkerSize= MSz*2;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = C{i} ;
    h(i).Color = 'k';%C(i);
    UL(i,:)=h(i).YData;
end
hold on
h=gscatter(double(catMethods),tmp.ul,tmp.blockContribution); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = '_';
    h(i).MarkerSize= MSz*2;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = C{i} ;
    h(i).Color = 'k';%C(i);
    for j=1:size(h(1).XData,2)
        line([h(i).XData(j),h(i).XData(j)],[UL(i,j),h(i).YData(j)],'Color','k')
    end
end
h=gscatter(double(catMethods),tmp.m,tmp.blockContribution); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = S{i};
    h(i).MarkerSize= MSz;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = 'k';%C(i) ;
    h(i).Color = 'k';%C(i);
end
t5.Title.String='DynOrder=="[2;0]"';
t5.Title.FontSize=26;
t5.YLim=[0 1];
t5.XLim=[0.5 5.5];
t5.XTick=1:5;
t5.XTickLabel=categories(catMethods);
t5.XTickLabelRotation=30;
t5.FontSize=18;
box on;

t7=nexttile(4);
tmp = (plotData(plotData.DynOrder=="[0;2]",:));%,:);%
catMethods=categorical(tmp.ModelType,varNames);
offset=[-0.15,-0.05,0.05,0.15];%[-0.1,0,0.1];
UL=nan(3,5);
h=gscatter(double(catMethods),tmp.ll,tmp.blockContribution); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = '_';
    h(i).MarkerSize= MSz*2;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = C{i} ;
    h(i).Color = 'k';%C(i);
    UL(i,:)=h(i).YData;
end
hold on
h=gscatter(double(catMethods),tmp.ul,tmp.blockContribution); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = '_';
    h(i).MarkerSize= MSz*2;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = C{i} ;
    h(i).Color = 'k';%C(i);
    for j=1:size(h(1).XData,2)
        line([h(i).XData(j),h(i).XData(j)],[UL(i,j),h(i).YData(j)],'Color','k')
    end
end
h=gscatter(double(catMethods),tmp.m,tmp.blockContribution); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = S{i};
    h(i).MarkerSize= MSz;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = 'k';%C(i) ;
    h(i).Color = 'k';%C(i);
end
t7.Title.String='DynOrder=="[0;2]"';
t7.Title.FontSize=26;
t7.YLim=[0 1];
t7.XLim=[0.5 5.5];
t7.XTick=1:5;
t7.XTickLabel=categories(catMethods);
t7.XTickLabelRotation=30;
t7.FontSize=18;
box on;

leg=legend(h,'Orientation', 'vertical');
leg.Layout.Tile = 'east';
leg.Interpreter = 'tex';
leg.Title.String="blockContribution";
box on;
T1.XLabel.String="Method";T1.XLabel.FontWeight="bold";T1.XLabel.FontSize=26;
T1.YLabel.String="R^2_P"+newline;T1.YLabel.FontWeight="bold";T1.YLabel.FontSize=26;