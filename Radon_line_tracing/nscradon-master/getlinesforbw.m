function bwlines = getlinesforbw(bw, params)
%GETLINESFORBW - detect lines on a binary image using Radon transform
%
% Input parameters:
% -----------------
% BW - binary image (it will be skeletonized in order to find lines)
% PARAMS - structure with main parameters
%  .THETA_STEP - resolution for theta values
%  .PEAK_THRESHOLD - threshold for detection of peaks on transform 
%  .LINE_FILLGAP - lines with this gaps will be merged into one 
%  .LINE_MINLENGTH - minimum length of line that will be detected 
%  .SPLIT_BW - split BW into segments and process segments individually or not
%  .SHOW_BW_PLOT - show control plots with peaks and detected lines
%  .ANGLETOL - maximum angle between two crossing lines to merge them into one
%
% Output:
% ------------------
% BWLINES - a structure array with detected lines (see RADONLINES for details)
%
% See also: RADONLINES, FINDPEAKS, SHOWFOUNDLINES
%
   
   % check parameters and set up default values
   if isfield(params, 'theta_step') 
      theta_step = params.theta_step;
   else
      theta_step = 0.1;
   end
         
   if isfield(params, 'peak_threshold')
      threshold = params.peak_threshold;
   else
      threshold = 0.3;
   end   
   
   if isfield(params, 'line_fillgap')
      line_fillgap = params.line_fillgap;
   else
      line_fillgap = 10;
   end   
   
   if isfield(params, 'line_minlength')
      line_minlength = params.line_minlength;
   else
      line_minlength = 10;
   end
   
   
   if isfield(params, 'split_bw')
      split_bw = params.split_bw;
   else
      split_bw = 0;
   end
   
   if isfield(params, 'show_bw_plots')
      show_bw_plots = params.show_bw_plots;
   else
      show_bw_plots = 0;
   end
   
   if isfield(params, 'angletol')
      angletol = params.angletol;
   else
      angletol = 0;
   end
   

   theta = 0:theta_step:180-theta_step;
   
   % split image into segments or process as one segment
   if split_bw == 1
      [L, n] = bwlabeln(bw);
   else
      L = bw;
      n = 1;
   end
      
   bwlines = [];
   for k = 1:n
      bw_seg = L == k;
               
      % skip small segments
      if sum(sum(bw_seg)) < line_minlength
         continue;
      end
               
      % detect lines using Radon transform
      [R, rho] = radon(bwmorph(bw_seg, 'skel', inf), theta);
      P = findpeaks(R, threshold);      
      lines = radonlines(bw_seg, theta, rho, P, line_fillgap, line_minlength, angletol);
      
      % accumulate detected lines and show control plots if needed   
      if numel(lines) > 0
         bwlines = [bwlines lines];
      end
      
      if show_bw_plots == 1 && k < 20
         show_plots(R, rho, theta, P, bw_seg, lines);
      end            
   end
end

function show_plots(R, rho, theta, P, bw, lines)
%SHOW_PLOTS shows conrol plots with peaks and lines
%
   figure
   subplot(4, 4, [1 5 9 13]);
   imagesc(theta, rho, R) 
   colormap(jet)
   ylabel('Rho');
   xlabel('Theta');
   title('Transform');
   
   subplot(4, 4, [2 6 10 14]);
   imagesc(theta, rho, R)
   colormap jet;
   ylabel('Rho');
   xlabel('Theta');
   colormap(jet);
   hold on
   for k = 1:size(P, 1)
      plot(theta(P(k, 2)), rho(P(k, 1)), 'og');
   end
   hold off
   title('Found peaks');

   subplot(4, 4, [3 4 7 8]);
   imshow(bw, []), colormap(gray);
   title('Binary image');
   
   subplot(4, 4, [11 12 15 16]);
   imshow(bw), colormap gray
   hold on
   for k = 1:numel(lines)
      p1 = lines(k).point1;
      p2 = lines(k).point2;
      line([p1(1) p2(1)], [p1(2) p2(2)], 'Color', 'red', 'LineWidth', 2);
   end
   hold off
   title(sprintf('Binary image (found %d lines)', numel(lines)));
end
   
