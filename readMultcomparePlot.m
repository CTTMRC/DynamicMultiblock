function T=readMultcomparePlot(f)
plotNames=f.Children(1).YTickLabel;
lineStack=f.Children(1).Children;
S1=split(plotNames,',');
[n,m]=size(S1);
ll=nan(n,1);
mm=nan(n,1);
ul=nan(n,1);
levelName=cell(m,1);
S2=struct;
for i=1:m
    tempSplit=split(S1(:,i),'=');
    levelName(i)=tempSplit(1,1);
    S2.(levelName{i})=string(tempSplit(:,2));
end
startIdx=[0:2:n*2-2]+1;
endIdx=2:2:n*2;
for i=1:n
    ll(i)=lineStack(startIdx(i)).XData(1);
    ul(i)=lineStack(startIdx(i)).XData(2);
    mm(i)=lineStack(endIdx(i)).XData;
end
T=table('Size',[n,m+3],'VariableTypes',[repmat("string",[m,1]);repmat("double",[3,1])]);
T.Properties.VariableNames=cellstr([string(fields(S2));"ll";"m";"ul"]); 
for i=fields(S2)'
    T.(string(i))=S2.(string(i));
end
T.ll=ll;
T.m=mm;
T.ul=ul;