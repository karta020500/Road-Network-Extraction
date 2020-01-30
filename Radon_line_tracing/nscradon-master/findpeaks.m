function P = findpeaks(R, threshold, gaussParam) 
%FINDPEAKS finds peaks on Hough or Radon transform result
%
   R = imfilter(R, fspecial('gaussian', 3, 2));   
   Rbw = R - imfilter(R, fspecial('gaussian', 80, 20)) > threshold * max(R(:));
   Rbw = imdilate(Rbw, strel('disk', 2));

   [L, n] = bwlabeln(Rbw);
   P = zeros(n, 2);
   for i = 1:n
      Ri = L == i;
      Ra = Ri .* R;
      Ramx = max(max(Ra));
      [r, c] = find(Ra == Ramx);
      idx = ceil(numel(r)/2);
      P(i, :) = [r(idx) c(idx)];
   end   
end
