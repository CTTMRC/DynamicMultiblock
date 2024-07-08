function fig=selectionHelperPlot(Q,R,varargin)
%Inputparser
defaultneedTilt=true;
checkInput=inputParser;
addRequired(checkInput, 'Q')
addRequired(checkInput, 'R',@(x)mustBeEqualSize(x,Q,1))
addOptional(checkInput, 'needTilt',defaultneedTilt,@islogical)
parse(checkInput,Q,R,varargin{:})
%init
Q=checkInput.Results.Q;
R=checkInput.Results.R;
needTilt=checkInput.Results.needTilt;
Q2yCI=cell(1,size(Q,3));
R2yCI=cell(1,size(R,3));
for k=1:size(Q,3)
    Q2yCI{k}={cumsum(bootci(15,@(x)mean(x),Q(:,:,k)),2)'};
    R2yCI{k}={cumsum(bootci(15,@(x)mean(x),R(:,:,k)),2)'};
end
if needTilt
    fig=figure('WindowState','maximized');
    g=tiledlayout(fig,size(Q,3),2);
    for i=1:size(Q,3)
        [M,I]=max((Q2yCI{i}{:,:}(:,1)'));
        t1=nexttile(g);
        scatter(t1,1:size(Q,2),cumsum(mean(Q(:,:,i))),'r+')
        hold on
        scatter(t1,1:size(Q,2),Q2yCI{i}{:,:},'k_')
        hold on
        patch(t1,[1:size(Q,2),size(Q,2):-1:1],[Q2yCI{i}{:,:}(:,1)',flipud(Q2yCI{i}{:,:}(:,2))'],...
            'r','FaceAlpha',0.1,'EdgeColor','r','EdgeAlpha',0.2)
        hold on
        scatter(t1,I,M,'gd','filled')
        hold on
        yline(t1,max(cumsum(mean(Q(:,:,i)))),'--',num2str(round(max(cumsum(mean(Q(:,:,i)))),2)),...
            'LabelVerticalAlignment','bottom','LineStyle','-.','Color','#818589','Alpha',0.5)
        YL=cat(2,Q2yCI{i});
        YL=cat(2,YL{:,:});
        t1YLim=[min(YL,[],'all')-0.1,min([max(YL,[],'all')+0.1,1])];
        if mod(i,2)
            t1.YLabel.String="% of Variance Explained"+newline;
        else
            t1.YLabel.String="% of Variance Explained";
        end
        t1.XLabel.String='Components';
        t1.Title.String='Q^2';
        t1.FontSize=10;
        box on
        
        t2=nexttile(g);
        scatter(t2,1:size(R,2),cumsum(mean(R(:,:,i))),'r+')
        hold on
        scatter(t2,1:size(R,2),R2yCI{i}{:,:},'k_')
        hold on
        patch(t2,[1:size(R,2),size(R,2):-1:1],[R2yCI{i}{:,:}(:,1)',flipud(R2yCI{i}{:,:}(:,2))'],'r','FaceAlpha',0.1,'EdgeColor','r','EdgeAlpha',0.2)
        YL=cat(2,R2yCI{i});
        YL=cat(2,YL{:,:});
        t2YLim=[min(YL,[],'all')-0.1,min([max(YL,[],'all')+0.1,1])];
        if mod(i,2)
            t2.YLabel.String="% of Variance Explained"+newline;
        else
            t2.YLabel.String="% of Variance Explained";
        end
        t2.XLabel.String='Components';
        t2.Title.String='R^2';
        t2.FontSize=10;
        box on
        t1.YLim=[min([t1YLim(1),t2YLim(1)]),max(([t1YLim(2),t2YLim(2)]))];
        t2.YLim=[min([t1YLim(1),t2YLim(1)]),max(([t1YLim(2),t2YLim(2)]))];
    end
elseif ~needTilt
    fig=figure('WindowState','maximized');
    g=tiledlayout(fig,2,size(Q,3));
    for i=1:size(Q,3)
        t1=nexttile(g,i);
        scatter(t1,1:size(Q,2),cumsum(mean(Q(:,:,i))),'r+')
        hold on
        scatter(t1,1:size(Q,2),Q2yCI{i},'k_')
        hold on
        patch(t1,[1:size(Q,2),size(Q,2):-1:1],[Q2yCI{i}(:,1)',flipud(Q2yCI{i}(:,2))'],'r','FaceAlpha',0.1,'EdgeColor','r','EdgeAlpha',0.2)
        box on
        t2=nexttile(g,i+size(Q,3));
        scatter(t2,1:size(R,2),cumsum(mean(R(:,:,i))),'r+')
        hold on
        scatter(t2,1:size(R,2),R2yCI{i},'k_')
        hold on
        patch(t2,[1:size(R,2),size(R,2):-1:1],[R2yCI{i}(:,1)',flipud(R2yCI{i}(:,2))'],'r','FaceAlpha',0.1,'EdgeColor','r','EdgeAlpha',0.2)
        box on
    end
else
    keyboard
end
%--------------------------------------------------------------------------
%   Helper functions
%--------------------------------------------------------------------------
% Custom validation function
function mustBeEqualSize(a,b,d)
% Test for equal size
if ~isequal(size(a,d),size(b,d))
    eid = 'Size:notEqual';
    msg = 'Size of first input must equal size of second input.';
    throwAsCaller(MException(eid,msg))
end