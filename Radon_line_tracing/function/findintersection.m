function [x y] = findintersection(L1,L2)

% findintersection is a function for finding the intersection point of two lines (vectors)
%   Inputs into this function are coordinates of two lines (vectors)
%   Line1 = [(x11,y11);(x12,y12)] and line2=[(x21,y21);(x22,y22)]
%   Example1 [xintersect yintersect] = findintersection([0.5 1;1 0],[0.5 0;1 1])
%   Example2 [xintersect yintersect]= findintersection([15 0;100 200],[50 0;50 200])
%   Example3 [xintersect yintersect]= findintersection([0 50;100 50],[50 0;50 100])
%   Outputs are x and y positions of the intersection point
%   Written by L M Vhengani
%   2011 September
%   Email: Lvhengani@gmail.com
%   Reference: http://mathworld.wolfram.com/Line-LineIntersection.html

x1=L1(1,1);y1=L1(1,2);x2=L1(2,1);y2=L1(2,2); % Components of Line 1
x3=L2(1,1);y3=L2(1,2);x4=L2(2,1);y4=L2(2,2); % Components of Line 2



x = det([det([x1 y1;x2 y2]), (x1-x2);det([x3 y3;x4 y4]), (x3-x4) ])/det([(x1-x2),(y1-y2) ;(x3-x4),(y3-y4)]);

y = det([det([x1 y1;x2 y2]), (y1-y2);det([x3 y3;x4 y4]), (y3-y4) ])/det([(x1-x2),(y1-y2) ;(x3-x4),(y3-y4)]);

x = fix(x);
y = fix(y);

end