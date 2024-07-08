function blockProgress(i,tot,init)
charRatio=tot/75;
if init==true
fprintf('\n');
else
fprintf(char(repmat(8,1,76)));
end
percent=ceil(i/charRatio);
fprintf([char(repmat(9632,1,percent)),char( repmat(9633,1,75-percent) ),'\n']);