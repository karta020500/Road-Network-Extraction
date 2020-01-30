function bwlines = radonlines(BW, thetas, rho, P, fillgap, minlength, angletol)
%RADONLINES detect lines on BW image using peaks found on Radon transform
%
% Input parameters:
% -----------------
% BW - binary image (it will be skeletonized in order to find lines)
% THETAS - a vector with values for Theta use to make transformation
% RHO - a vector with values for Rho use to make transformation
% FILLGAP - lines with this or smaller gaps will be merged into one 
% MINLENGTH - minimum length of line that will be detected 
% ANGLETOL - if two lines cross each other with angle smaller than this, merge them
%
% Output:
% -------
% BWLINES - a structure array with detected lines, inlucindg
%  .POINT1 - vector with coordinates of point 1 of a line
%  .POINT2 - vector with coordinates of point 2 of a line
%  .LENGTH - length of the line (in pixels)
%
% See also: FINDPEAKS, GETLINESFORBW, SHOWFOUNDLINES
%
   
   if nargin < 7
      angletol = 4;
   end
   
   if nargin < 6
      minlength = 10;
   end
   
   if nargin < 5
      fillgap = 10;
   end
   
   % square the distance related values to avoid use of sqrt()
   fillgap_sq = fillgap^2;
   minlength_sq = minlength^2;
   
   % find the center and make coordinates transformation
   [h, w] = size(BW);
   x0 = ceil(w/2);
   y0 = ceil(h/2);

   [y, x] = find(BW);
   nonzeropix = [x, y];
   
   x = nonzeropix(:,1);
   y = nonzeropix(:,2);
   xx = x - x0;
   yy = y - y0;
   
   % loop for every peak found
   numlines = 0;
   for k = 1:size(P, 1)
      
      % find non-zeros BW pixels that are located on the direction
      % for given theta and rho
      theta = -thetas(P(k, 2));   
      theta_c = theta * pi / 180;
      rho_xy = xx * cos(theta_c) + yy * sin(theta_c);
      idx = round(rho_xy) == rho(P(k, 1));      
      r = y(idx) + 1; 
      c = x(idx) + 1;      
      
      % resort the pixels according to the direction
      [r, c] = reSortHoughPixels(r, c);      
   
      if isempty(r) 
         continue 
      end
      
      % compute distance^2 between the point pairs
      xy = [c r]; % x,y pairs in coordinate system with the origin at (1,1)
      diff_xy_sq = diff(xy,1,1).^2;
      dist_sq = sum(diff_xy_sq,2);
      
      % find the gaps between points larger than the FILLGAP
      fillgap_idx = find(dist_sq > fillgap_sq);
      idx = [0; fillgap_idx; size(xy,1)];
      
      % loop for the pixels 
      for p = 1:length(idx) - 1
         p1 = xy(idx(p) + 1, :); % offset by 1 to convert to 1 based index
         p2 = xy(idx(p + 1), :); % set the end (don't offset by one this time)

         % if line length is larger than MINLENGTH, save start and end
         % coordinates
         linelength_sq = sum((p2 - p1).^2);
         if linelength_sq >= minlength_sq
            numlines = numlines + 1;
            bwlines(numlines).point1 = p1;
            bwlines(numlines).point2 = p2;
         end
      end
   end   
   
   if numlines == 0
      bwlines = [];
      return;
   end   
      
   % find lines which intersect each other with angle smaller than
   % ANGLETOL and merge them into one
   if angletol > 0
      n = 1;
      while n <= numel(bwlines)
         l1 = bwlines(n);
         q = n + 1;
         while q <= numel(bwlines)
            l2 = bwlines(q);
            % x-y coordinates of the end points
            x = [l1.point1(1), l1.point2(1), l2.point1(1), l2.point2(1)];
            y = [l1.point1(2), l1.point2(2), l2.point1(2), l2.point2(2)];
            
            % crossing point
            den = (x(1) - x(2))*(y(3) - y(4)) - (y(1) - y(2))*(x(3) - x(4));
            px = (x(1)*y(2) - y(1)*x(2))*(x(3) - x(4)) - (x(1) - x(2))*(x(3)*y(4) - y(3)*x(4));
            py = (x(1)*y(2) - y(1)*x(2))*(y(3) - y(4)) - (y(1) - y(2))*(x(3)*y(4) - y(3)*x(4));
            px = px / den;
            py = py / den;
            
            % slope of each line
            m1 = (y(2) - y(1)) / (x(2) - x(1));
            m2 = (y(4) - y(3)) / (x(4) - x(3));
            
            % angle between lines
            alpha = atand( (m1 - m2) / (1 + m1 * m2) );

            % sort points
            [xx, ind] = sort(x);
            yy = y(ind);

            % get pairwise distance between line points
            xy = [x', y'];
            pairs = nchoosek(1:4, 2);
            p1 = xy(pairs(:,1),:);
            p2 = xy(pairs(:,2),:);
            d = sqrt(sum((p1-p2).^2,2));
            
            debug = 0;
            if debug
               disp('------')
               disp([x, px])
               disp([y, py])
               disp(alpha)
               disp(d)
               figure
               hold on
               plot(x(1:2), y(1:2), '.-r')
               plot(x(3:4), y(3:4), '.-b')
               scatter(px, py, 'xk')
               hold off
               title([alpha, min(d)])
            end
            
            % if
            % 1. angle between lines is smaller than ANGLETOL
            % 2. lines are close to each other or intersection point is within the area spanned 
            %    by both lines
            % merge the lines by taking average of end points            
            if abs(alpha) < angletol && ...
               (min(d) <= fillgap ||...
                  (px >= min(xx) - fillgap && px <= max(xx) + fillgap && ...
                  py >= min(yy) - fillgap && py <= max(yy) + fillgap))...
                
               
            
               bwlines(n).point1 = [xx(1); yy(1)];
               bwlines(n).point2 = [xx(4); yy(4)];
               bwlines(q) = [];
               if debug
                  disp('merged!')
                  hold on
                  plot([bwlines(n).point1(1) bwlines(n).point2(1)], ...
                     [bwlines(n).point1(2) bwlines(n).point2(2)], '--k')
                  hold off
                  title('merged!')
               end
            else               
               q = q + 1;               
            end
         end
         bwlines(n).length = sqrt(sum((bwlines(n).point2 - bwlines(n).point1).^2));
         n = n + 1;
      end   
   end 
end

function [r_new, c_new] = reSortHoughPixels(r, c)
% make sure that r an c are in the order along the line segment

   if isempty(r)
      r_new = r;
      c_new = c;
      return;
   end

   r_range = max(r) - min(r);
   c_range = max(c) - min(c);
   
   if r_range > c_range
      % Sort first on r, then on c
      sorting_order = [1 2];
   else
      % Sort first on c, then on r
      sorting_order = [2 1];
   end

   [rc_new] = sortrows([r c], sorting_order);
   r_new = rc_new(:,1);
   c_new = rc_new(:,2);
end
