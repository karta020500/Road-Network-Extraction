function showfoundlines(bw, bwlines)
% SHOWFOUNDLINES shows an image and lines detected on it
%
   imshow(bw)
   hold on
   for l = 1:numel(bwlines)
      p1 = bwlines(l).point1;
      p2 = bwlines(l).point2;
      line([p1(1) p2(1)], [p1(2) p2(2)], 'Color', 'r', 'LineWidth', 1);
   end  
   hold off
end