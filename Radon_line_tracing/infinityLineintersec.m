function point = infinityLineintersec(x1,y1,x2,y2)

% x1  = [1345 1345];
% y1  = [1003 1102];
% %line2
% x2 = [1252 1445];
% y2 = [1134 1137];
%fit linear polynomial
p1 = polyfit(x1,y1,1);
p2 = polyfit(x2,y2,1);
%calculate intersection
x_intersect = fzero(@(x) polyval(p1-p2,x),3);
y_intersect = polyval(p1,x_intersect);
point = [x_intersect,y_intersect];