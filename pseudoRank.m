function pr=pseudoRank(x)
oldFolder=cd([toolboxdir('stats'),'\stats']);
[~,~,~,~,PR,~]=pca(x);
pr=find(cumsum(PR)>90,1,'first');
cd(oldFolder);
end