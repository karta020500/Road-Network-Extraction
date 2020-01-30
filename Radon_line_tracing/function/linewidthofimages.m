function [f, seg_list_new] = linewidthofimages(imagepath)
%bundleImage=imread(imagepath);
[bundleImage,R] = geotiffread(imagepath);
R2 = worldfileread('Road_calss.tfw', 'planar', size(bundleImage));

bundleImage(bundleImage > 3) = 0;
bundleImage(bundleImage == 3) = 1;

% % convert to greyscale if necessary
% if  size(bundleImage,3)>1
%     bundleImage = rgb2gray(bundleImage);
% end;

  text_image = zeros(500);

%H
% text_image(250:260,:) = 1;

%V
%text_image(:,250:260) = 1;

%D
% text_image = eye(500);
% text_image = imdilate(text_image,ones(3,3));

%test
% bundleImage = text_image;
% 
[m,n] = size(bundleImage);

%crop
% if ~exist('cropRect')
%     cropRect=[10,10,size(bundleImage,2)-20,size(bundleImage,1)-20];
% end;
% f=figure;
%
% imshow(bundleImage,'Border','tight');
% set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
% h=imrect(gca,cropRect)
% cropRect = wait(h);
% [bundleImage]=imcrop(bundleImage,cropRect);
% close(f)

myfigure = figure('Name',imagepath,'NumberTitle', 'off');
bundleForPeakFinding = bundleImage;

%subplot(2,2,1); imshow(bundleImage);



bundleImage=double(bundleImage);
subplot(2,2,1); mapshow(bundleImage,R2);
imagesc(bundleImage);
axis image off

%pre-processing------------------------------------------------------------

%opening
se = strel('square',2);
opened = imopen(bundleImage,se);

%closeing
se2 = strel('square',5);
closed = imclose(opened,se2);

%connectivity_area
Pixel_connectivity = bwareaopen(closed, 50);


%hole_fill
filled = imfill(Pixel_connectivity,'holes');
holes = filled & ~Pixel_connectivity;
bigholes = bwareaopen(holes, 200);
smallholes = holes & ~bigholes;
bundleImage = Pixel_connectivity | smallholes;
% 
skel = bwmorph(filled,'skel',Inf);

skel = bwmorph(skel, 'spur', 9);
% 
skel = imdilate(skel,ones(5,5));
%-------------------------------------------------------------

% iptsetpref('ImshowAxesVisible','on')
% theta = 0:179;
% [R,xp] = radon(bundleForPeakFinding,theta);

% subplot(2,2,2); imshow(R,[],'Xdata',theta,'Ydata',xp,'InitialMagnification','fit')
% xlabel('\theta (degrees)')
% ylabel('x''')
% colormap(hot), colorbar
% axis normal
% iptsetpref('ImshowAxesVisible','off')

% peaks=peakfit2d(R);
% peaks=round(peaks);
% 
% subplot(2,2,3); imagesc(R); hold on
% subplot(2,2,3); plot(peaks(2),peaks(1),'r+')
% subplot(2,2,4); plot(xp,R(:,peaks(2)));
% title('R_{0^o} (x\prime)')
% 
% [R,xp] = radon(bundleImage,peaks(2));
% width1=fwhm(xp,R);
% 
% stringout=sprintf('Width:%0.2f',width1);
% disp(stringout)
% 
% f= width1;
% 
% 


% theta = [0:180]';
theta_step = 0.5;
theta = [0:theta_step:180-theta_step]';
[R, rho] = radon(bundleImage, theta) ;
 
%R = medfilt2(R,[3,3]);
R = imfilter(R, fspecial('gaussian', 3, 2));

subplot(2,2,2); imshow(R,[],'Xdata',theta,'Ydata',rho,'InitialMagnification','fit')
xlabel('\theta'), ylabel('\rho');
colormap(hot), colorbar
axis on, axis normal
iptsetpref('ImshowAxesVisible','off')

subplot(2,2,3) ;
imshow(R, [], 'XData', theta, 'YData', rho, 'InitialMagnification', 'fit');
xlabel('\theta'), ylabel('\rho');
axis on, axis normal, hold on;
%-------------------------------------------------------------
% Detect the peaks of transform output  
% Threshold value for peak detection
threshold_val = ceil(0.40*max(R(:)));

% initialize the loop variables
done = false;
h = R;
hnew = h;
nhood = size(R)/10;
nhood = max(2*ceil(nhood/2) + 1, 1);
nhood_center = (nhood-1)/2;
P = [];
numpeaks = 100;

while ~done
  [dummy max_idx] = max(hnew(:)); %#ok
  [p, q] = ind2sub(size(hnew), max_idx);
  
  p = p(1); q = q(1);
  if hnew(p, q) >= threshold_val
    P = [P; [p q]]; %#ok<AGROW> % add the peak to the list
    
    % Suppress this maximum and its close neighbors.
    p1 = p - nhood_center(1); p2 = p + nhood_center(1);
    q1 = q - nhood_center(2); q2 = q + nhood_center(2);
    % Throw away neighbor coordinates that are out of bounds in
    % the rho direction.
    [qq, pp] = meshgrid(q1:q2, max(p1,1):min(p2,size(h,1)));
    pp = pp(:); qq = qq(:);
    
    % For coordinates that are out of bounds in the theta
    % direction, we want to consider that H is antisymmetric
    % along the rho axis for theta = +/- 90 degrees.
    theta_too_low = find(qq < 1);
    qq(theta_too_low) = size(h, 2) + qq(theta_too_low);
    pp(theta_too_low) = size(h, 1) - pp(theta_too_low) + 1;
    theta_too_high = find(qq > size(h, 2));
    qq(theta_too_high) = qq(theta_too_high) - size(h, 2);
    pp(theta_too_high) = size(h, 1) - pp(theta_too_high) + 1;
    
    % Convert to linear indices to zero out all the values.
    hnew(sub2ind(size(hnew), pp, qq)) = 0;
    
    done = size(P,1) == numpeaks;
  else
    done = true;
  end
end

P = unique(P,'rows');

% % Maximum nof peaks to identify in the image
% 
% max_nofpeaks = 1;
% max_indexes = find(R(:)>threshold_val) ;
% 
% max_values = R(max_indexes) ;
% [sorted_max, temp_indexes] = sort(max_values, 'descend') ;
% sorted_indexes = max_indexes(temp_indexes) ;
%     
% 
% % Get the first highest peaks for the sorted array
% if (length(sorted_max) <= max_nofpeaks)
%     peak_values = sorted_max(1:end) ; 
%     peak_indexes = sorted_indexes(1:end) ;
% else
%     peak_values = sorted_max(1:max_nofpeaks) ;
%     peak_indexes = sorted_indexes(1:max_nofpeaks) ;
% end
% [y, x]  = ind2sub(size(R), peak_indexes ) ;
% peaks = [rho(y) theta(x)] ;

% Fast_P = FastPeakFind(R);
% 
% Px_index = Fast_P(1:2:end);
% Py_index = Fast_P(2:2:end);
% 
% P(:,1) = Fast_P(2:2:end);
% P(:,2) = Fast_P(1:2:end);
% 
% peaks = [rho(Py_index) theta(Px_index)];

% Peak profile analysis

[p_m, p_n] = size(P);
[r_m, r_n] = size(R);
range = 5;
semi_peak_radon = 0.9;
Width_seg={};
count=1;

%linear interpolation
for i = 1:p_m
    query = R(P(i,1),P(i,2))*semi_peak_radon;
    
    if (P(i,1) > range && P(i,1) < r_m-range) && (P(i,2) > range && P(i,2) < r_n-range)
       
       
       %offset_p conputation
       p_range = R(P(i,1)-range:P(i,1)+range,P(i,2));
       [~,p_range_I]=sort(abs(p_range-query));
      
%        p_query = (p_range(p_range_I(1)) + p_range(p_range_I(2)))/2;
       p_width = p_range(p_range_I(1)) - p_range(p_range_I(2));
%        [~,p_query_I] = min(abs(p_range - p_query));
%        
%        p0 = p_range(p_query_I); 

       
       %theta_the conputation(interpolation)
%        the_range_the = R(P(i,1),P(i,2)-range:P(i,2)+range);
%        the_range_p =   theta(P(i,2)-range:P(i,2)+range)';
%        interp_intervel = (min(the_range_p):.1:max(the_range_p));
%        
%        range_interp = interp1q(the_range_p',the_range_the',interp_intervel');
%        [~,range_I]=sort(abs(range_interp-query));
%        
%        the_query = (range_interp(range_I(1)) + range_interp(range_I(2)))/2;
%        [~,the_query_I] = min(abs(range_interp - the_query));
%        
%        the0 = range_interp(the_query_I);

       %Replace original values
       %P(i,1) = P(i,1) + ((range+1) - p_query_I);
       %P(i,2) = P(i,2) + ((range*10+1) - the_query_I)/10;
       Width_seg{count}(1,1) = i;
       Width_seg{count}(1,2) = abs(p_width);
       count = count+1;
       
    end
end

peaks = [rho(P(:,1)) theta(P(:,2))];
plot(peaks(:,2), peaks(:,1), 's', 'color','white');
title('Radon Transform & Peaks') ;

subplot(2,2,4); plot(rho,R(:,P(1,2)));

[R_w,rho_w] = radon(bundleImage,P(2));
width1=fwhm(rho_w,R_w);

stringout=sprintf('Width:%0.2f',width1);
disp(stringout)

f= width1;
%-------------------------------------------------------------

% Detected lines on the image
figure()
imshow(bundleImage), title('Detected lines'), hold on

x_center = floor(size(bundleImage, 2)/2) ;
y_center = floor(size(bundleImage, 1)/2) ;

x1 = [];
y1 = [];
x2 = [];
y2 = [];

for p=1:p_m

    x_1 = [-x_center, x_center];
    y_1 = [0, 0] ;

    % Shift at first
    x_1_shifted = x_1 ;
    y_1_shifted = [y_1(1)-peaks(p,1), y_1(2)-peaks(p,1)];

    % Rotate 
    peaks(p,2) = 90 - peaks(p,2);
    t=peaks(p,2)*pi/180;
    rotation_mat = [ cos(t) -sin(t) ; sin(t) cos(t) ];
    x_y_rotated = rotation_mat*[x_1_shifted; y_1_shifted];
    x_rotated = x_y_rotated(1,:) ;
    y_rotated = x_y_rotated(2,:) ;
    
    x1(p) = x_rotated(1)+x_center;
    y1(p) = y_rotated(1)+y_center;
    x2(p) = x_rotated(2)+x_center;
    y2(p) = y_rotated(2)+y_center;
    
    plot([x1(p) x2(p)], [y1(p) y2(p)], 'r', 'linewidth', 2);
end

hold off;

%-------------------------------------------------------------

% Create a line from point 1 to point 2
% line_m = zeros(m,n);
seg_list={};
% for l=1:p_m
%     spacing = 1;
%     line_s = zeros(m,n);
%     numSamples = ceil(sqrt((x2(l)-x1(l))^2+(y2(l)-y1(l))^2) / spacing);
%     x = linspace(x1(l), x2(l), numSamples);
%     y = linspace(y1(l), y2(l), numSamples);
%     xy = round([x',y']);
%     dxy = abs(diff(xy, 1));
%     duplicateRows = [0; sum(dxy, 2) == 0];
%     % Now for the answer:
%     finalxy = xy(~duplicateRows,:);
%     finalx = finalxy(:, 1);
%     finaly = finalxy(:, 2);
%     % Plot the points
%     hold on;
%     %plot(finalx, finaly, 'y*');
%     
%     
%     for p=1:length(finalxy)
%         if  finalxy(p,1) > 0 && finalxy(p,1) < n && finalxy(p,2) > 0 && finalxy(p,2) < m
%         line_m(finaly(p),finalx(p)) = bundleImage(finaly(p),finalx(p));
%         line_s(finaly(p),finalx(p)) = bundleImage(finaly(p),finalx(p));
%         end
%     end
%     
%     line_seg = find(line_s);
%     [I,J] = ind2sub(size(line_s),max(line_seg));
%     [K,F] = ind2sub(size(line_s),min(line_seg));
%    
%     seg_list{l}(1,1) = I; seg_list{l}(1,2) = J;
%     seg_list{l}(2,1) = K; seg_list{l}(2,2) = F;
%     hold off;
% end

lines = radonlines(bundleImage,theta,rho,P,30,60,0);

for k = 1:numel(lines)
  seg_list{k}(1,1) = lines(k).point1(2); seg_list{k}(1,2) = lines(k).point1(1);
  seg_list{k}(2,1) = lines(k).point2(2); seg_list{k}(2,2) = lines(k).point2(1);
  %line([p1(1) p2(1)], [p1(2) p2(2)], 'Color', 'red', 'LineWidth', 2);
end


seg_index = [];
inter_x = [];
inter_y = [];
    for i = 2:length(seg_list)
        for j = 1:length(seg_list)
            if i == j && j < length(seg_list)
               j=j+1;
            end
            dt1=det([1,1,1;seg_list{i}(1,2),seg_list{i}(2,2),seg_list{j}(1,2);seg_list{i}(1,1),seg_list{i}(2,1),seg_list{j}(1,1)])*det([1,1,1;seg_list{i}(1,2),seg_list{i}(2,2),seg_list{j}(2,2);seg_list{i}(1,1),seg_list{i}(2,1),seg_list{j}(2,1)]);
            dt2=det([1,1,1;seg_list{i}(1,2),seg_list{j}(1,2),seg_list{j}(2,2);seg_list{i}(1,1),seg_list{j}(1,1),seg_list{j}(2,1)])*det([1,1,1;seg_list{i}(2,2),seg_list{j}(1,2),seg_list{j}(2,2);seg_list{i}(2,1),seg_list{j}(1,1),seg_list{j}(2,1)]);

            if(dt1<=0 && dt2<=0)
              intrsct = 1;         %If lines intesect
            else
              intrsct = 0;
            end
            
            if intrsct == 1
                diff_ang = (atan((seg_list{i}(2,1)-seg_list{i}(1,1))/(seg_list{i}(2,2)-seg_list{i}(1,2))) - atan((seg_list{j}(2,1)-seg_list{j}(1,1))/(seg_list{j}(2,2)-seg_list{j}(1,2)))) * 180/pi;    
                if abs(diff_ang) < 30
                [inter_x,inter_y] = polyxpoly(seg_list{i}(:,2),seg_list{i}(:,1),seg_list{j}(:,2),seg_list{j}(:,1));
                %seg_index = [seg_index;j];
                dist_i=bsxfun(@hypot,seg_list{i}(:,1)-inter_y,seg_list{i}(:,2)-inter_x);
                dist_j=bsxfun(@hypot,seg_list{j}(:,1)-inter_y,seg_list{j}(:,2)-inter_x);
                out_i = seg_list{i}(dist_i==min(dist_i),:);
                out_j = seg_list{j}(dist_j==min(dist_j),:);
                seg_list{i}(find(seg_list{i}==out_i)) = [inter_y,inter_x];
                seg_list{j}(find(seg_list{j}==out_j)) = [inter_y,inter_x];
%                 if true(seg_list_new{i}(1,:)==out_i)
%                    seg_list_new{i}(1,:) = [inter_y,inter_x];
%                 else
%                    seg_list_new{i}(2,:) = [inter_y,inter_x];
%                 end
% 
%                 if true(seg_list_new{j}(1,:)==out_j)
%                    seg_list_new{j}(1,:) = [inter_y,inter_x];
%                 else
%                    seg_list_new{j}(2,:) = [inter_y,inter_x];
%                 end
                end
            end
        end
    end    
    %seg_list(seg_index) =[];
    
%linkage
add_seg = {};
count = 1;

add_seg{count} = [1153, 731; 569,741];%street tree
count = count+1;


for i = 1:length(seg_list)
    for j = 1:length(seg_list)
         Distance = DistBetween2Segment(seg_list{i}(1, :), seg_list{i}(2, :), seg_list{j}(1, :), seg_list{j}(2, :));
%Distance between middle point of two line segments
         if abs(Distance)> 0 && abs(Distance) < 200
            diff_ang_seg = (atan((seg_list{i}(2,1)-seg_list{i}(1,1))/(seg_list{i}(2,2)-seg_list{i}(1,2))) - atan((seg_list{j}(2,1)-seg_list{j}(1,1))/(seg_list{j}(2,2)-seg_list{j}(1,2)))) * 180/pi;
            if abs(diff_ang_seg) >= 160 || abs(diff_ang_seg) <=20
               dist=bsxfun(@hypot,seg_list{j}(:,1)-seg_list{i}(1,1),seg_list{j}(:,2)-seg_list{i}(1,2));
               out = seg_list{j}(dist==min(dist),:);
               dist=bsxfun(@hypot,seg_list{i}(:,1)-out(1,1),seg_list{i}(:,2)-out(1,2));
               final = seg_list{i}(find(dist==min(dist), 1),:);
               diff_ang_conn_li = (atan((seg_list{i}(2,1)-seg_list{i}(1,1))/(seg_list{i}(2,2)-seg_list{i}(1,2))) - atan((out(1,1)-final(1,1))/(out(1,2)-final(1,2)))) * 180/pi;
%Lateral offset and direction difference of line segment
               if abs(diff_ang_conn_li) >= 160 || abs(diff_ang_conn_li) <=20
                  d = pdist([out(1,:); final(1,:)],'euclidean');
                  if abs(d)> 0 && abs(d) < 200
                     add_seg{count}(1,:) = out;
                     add_seg{count}(2,:) = final;
                     count = count+1;
                  end
               end
            end
        end
    end
end

 for i = 1:length(add_seg)  
      seg_list{length(seg_list)+1} = add_seg{i};
 end
     line_m = zeros(m,n);
    for i=1:length(seg_list)
        spacing = 0.5;
        numSamples = ceil(sqrt((seg_list{i}(2,2)-seg_list{i}(1,2))^2+(seg_list{i}(2,1)-seg_list{i}(1,1))^2) / spacing);
        x = linspace(seg_list{i}(1,2), seg_list{i}(2,2), numSamples);
        y = linspace(seg_list{i}(1,1), seg_list{i}(2,1), numSamples);
        xy = round([x',y']);
        dxy = abs(diff(xy, 1));
        duplicateRows = [0; sum(dxy, 2) == 0];
        % Now for the answer:
        finalxy = xy(~duplicateRows,:);
        finalx = finalxy(:, 1);
        finaly = finalxy(:, 2);
        % Plot the points

        for p=1:length(finalxy)
            if  finalxy(p,1) > 0 && finalxy(p,1) < n && finalxy(p,2) > 0 && finalxy(p,2) < m
            line_m(finaly(p),finalx(p)) = 1;
            end
        end
    end
    
line_m = filledgegaps(line_m, 3);
[edgelist, labelededgeim] = edgelink(line_m, 1);
seg_list_new = lineseg(edgelist, 10);
%drawedgelist(seg_list_new, size(line_m), 2, 'rand', 3); axis off

% %Junction connection-------------------------------------------------------
% r = 100;
% add_seg = {};
% count = 1;
% for i = 1:length(seg_list_new)
%      for j = 1:length(seg_list_new)
%           if i == j && j < length(seg_list_new)
%                j=j+1;
%           end
%           [po_x,po_y] = polyxpoly(seg_list_new{i}(:,2),seg_list_new{i}(:,1),seg_list_new{j}(:,2),seg_list_new{j}(:,1));
%           if isempty(po_x)
%              P = [seg_list_new{i}(1,2), seg_list_new{i}(1,1)];
%              line_seg = [seg_list_new{j}(1:end,2) seg_list_new{j}(1:end,1)];
%              polyout = polybuffer(P,'points',r);
%              [in,out] = intersect(polyout,line_seg);
%              if ~isempty(in)
%                 [interx,intery] = findintersection([seg_list_new{i}(end-1,2),seg_list_new{i}(end-1,1);seg_list_new{i}(end,2),seg_list_new{i}(end,1)],[seg_list_new{j}(end-1,2),seg_list_new{j}(end-1,1);seg_list_new{j}(end,2),seg_list_new{j}(end,1)]);
%                 point = [interx,intery];
%                  if point(1,1)<=1489 && point(1,2)<=1171
%                     segment_i = [seg_list_new{i}(1:end,2)  seg_list_new{i}(1:end,1)];
%                     segment_j = [seg_list_new{j}(1:end,2)  seg_list_new{j}(1:end,1)];
%                     polyout_2 = polybuffer(point,'points',1);
%                     [checkPt_i,out] = intersect(polyout_2,segment_i);
%                     [checkPt_j,out] = intersect(polyout_2,segment_j);
%                  if ~isempty(checkPt_i)
%                     dist=bsxfun(@hypot,seg_list_new{j}(:,2)-point(1,1),seg_list_new{j}(:,1)-point(1,2));
%                     if dist(1) == dist(2)
%                        dist = unique(dist);
%                     end
%                     short_end = seg_list_new{j}(dist==min(dist),:);
%                     %seg_list_new{j}(dist==min(dist),:) = point;
%                     d = pdist([short_end(1,2), short_end(1,1); point(1,:)],'euclidean');
%                     if abs(d)> 0 && abs(d) < r
%                        add_seg{count}(1,:) = short_end;
%                        add_seg{count}(2,:) = [point(1,2),point(1,1)];
%                        count = count+1;
%                     end
%                  elseif ~isempty(checkPt_j)
%                         dist=bsxfun(@hypot,seg_list_new{i}(:,2)-point(1,1),seg_list_new{i}(:,1)-point(1,2));
%                         if dist(1) == dist(2)
%                            dist = unique(dist);
%                         end
%                         short_end = seg_list_new{i}(dist==min(dist),:);
%                         %seg_list_new{i}(dist==min(dist),:) = point;
%                         d = pdist([short_end(1,2), short_end(1,1); point(1,:)],'euclidean');
%                         if abs(d)> 0 && abs(d) < r
%                            add_seg{count}(1,:) = short_end;
%                            add_seg{count}(2,:) = [point(1,2),point(1,1)];
%                            count = count+1;
%                         end
% %                  else
% %                        dist=bsxfun(@hypot,seg_list_new{i}(:,2)-point(1,1),seg_list_new{i}(:,1)-point(1,2));
% %                        if dist(1) == dist(2)
% %                           dist = unique(dist);
% %                        end
% %                        short_end = seg_list_new{i}(dist==min(dist),:);
% %                        %seg_list_new{i}(dist==min(dist),:) = point;
% %                        d = pdist([short_end(1,2), short_end(1,1); point(1,:)],'euclidean');
% %                        if abs(d)> 0 && abs(d) < r
% %                           add_seg{count}(1,:) = short_end;
% %                           add_seg{count}(2,:) = [point(1,2),point(1,1)];
% %                           count = count+1;
% %                        end
% %                        dist=bsxfun(@hypot,seg_list_new{j}(:,2)-point(1,1),seg_list_new{j}(:,1)-point(1,2));
% %                        if dist(1) == dist(2)
% %                           dist = unique(dist);
% %                        end
% %                        short_end = seg_list_new{j}(dist==min(dist),:);
% %                        %seg_list_new{j}(dist==min(dist),:) = point;
% %                        d = pdist([short_end(1,2), short_end(1,1); point(1,:)],'euclidean');
% %                        if abs(d)> 0 && abs(d) < r
% %                           add_seg{count}(1,:) = short_end;
% %                           add_seg{count}(2,:) = [point(1,2),point(1,1)];
% %                           count = count+1;
% %                        end
%                  end
%                  end
%              end
%              P = [seg_list_new{i}(end,2), seg_list_new{i}(end,1)];
%              line_seg = [seg_list_new{j}(1:end,2) seg_list_new{j}(1:end,1)];
%              polyout = polybuffer(P,'points',r);
%              [in,out] = intersect(polyout,line_seg);
%              if ~isempty(in)
%                 [interx,intery] = findintersection([seg_list_new{i}(end-1,2),seg_list_new{i}(end-1,1);seg_list_new{i}(end,2),seg_list_new{i}(end,1)],[seg_list_new{j}(end-1,2),seg_list_new{j}(end-1,1);seg_list_new{j}(end,2),seg_list_new{j}(end,1)]);
%                 point = [interx,intery];
%                  if point(1,1)<=1489 && point(1,2)<=1171
%                    segment_i = [seg_list_new{i}(1:end,2)  seg_list_new{i}(1:end,1)];
%                    segment_j = [seg_list_new{j}(1:end,2)  seg_list_new{j}(1:end,1)];
%   %                  [checkPt_i, onEnd] = checkPointOnSegment(segment_i,point,0);
%   %                  [checkPt_j, onEnd] = checkPointOnSegment(segment_j,point,0);
%                    polyout_2 = polybuffer(point,'points',1);
%                    [checkPt_i,out] = intersect(polyout_2,segment_i);
%                    [checkPt_j,out] = intersect(polyout_2,segment_j);
%                  if checkPt_i
%                     dist=bsxfun(@hypot,seg_list_new{j}(:,2)-point(1,1),seg_list_new{j}(:,1)-point(1,2));
%                     if dist(1) == dist(2)
%                        dist = unique(dist);
%                     end
%                     short_end = seg_list_new{j}(dist==min(dist),:);
%                     %seg_list_new{j}(dist==min(dist),:) = point;
%                     d = pdist([short_end(1,2), short_end(1,1); point(1,:)],'euclidean');
%                     if abs(d)> 0 && abs(d) < r
%                        add_seg{count}(1,:) = short_end;
%                        add_seg{count}(2,:) = [point(1,2),point(1,1)];
%                        count = count+1;
%                     end
%                  elseif checkPt_j
%                         dist=bsxfun(@hypot,seg_list_new{i}(:,2)-point(1,1),seg_list_new{i}(:,1)-point(1,2));
%                         if dist(1) == dist(2)
%                            dist = unique(dist);
%                         end
%                         short_end = seg_list_new{i}(dist==min(dist),:);
%                         %seg_list_new{i}(dist==min(dist),:) = point;
%                         d = pdist([short_end(1,2), short_end(1,1); point(1,:)],'euclidean');
%                         if abs(d)> 0 && abs(d) < r
%                            add_seg{count}(1,:) = short_end;
%                            add_seg{count}(2,:) = [point(1,2),point(1,1)];
%                            count = count+1;
%                         end
% %                  else
% %                        dist=bsxfun(@hypot,seg_list_new{i}(:,2)-point(1,1),seg_list_new{i}(:,1)-point(1,2));
% %                        if dist(1) == dist(2)
% %                           dist = unique(dist);
% %                        end
% %                        short_end = seg_list_new{i}(dist==min(dist),:);
% %                        %seg_list_new{i}(dist==min(dist),:) = point;
% %                        d = pdist([short_end(1,2), short_end(1,1); point(1,:)],'euclidean');
% %                        if abs(d)> 0 && abs(d) < r
% %                           add_seg{count}(1,:) = short_end;
% %                           add_seg{count}(2,:) = [point(1,2),point(1,1)];
% %                           count = count+1;
% %                        end
% %                        dist=bsxfun(@hypot,seg_list_new{j}(:,2)-point(1,1),seg_list_new{j}(:,1)-point(1,2));
% %                        if dist(1) == dist(2)
% %                           dist = unique(dist);
% %                        end
% %                        short_end = seg_list_new{j}(dist==min(dist),:);
% %                        %seg_list_new{j}(dist==min(dist),:) = point;
% %                        d = pdist([short_end(1,2), short_end(1,1); point(1,:)],'euclidean');
% %                        if abs(d)> 0 && abs(d) < r
% %                           add_seg{count}(1,:) = short_end;
% %                           add_seg{count}(2,:) = [point(1,2),point(1,1)];
% %                           count = count+1;
% %                        end
%                  end
%                  end
%              end
%           end 
%      end
% end

add_seg{count} = [1069,185;1135,186];

for i = 1:length(add_seg)  
     seg_list_new{length(seg_list_new)+1} = add_seg{i};
end

remove = [177,131;218,133];%short
index = cellfun(@(xx)isequal(xx,remove),seg_list_new);
[row,col] = find(index);
if ~isempty(col)
   for i = 1:length(col)
     seg_list_new(col(i)) = [];
   end
end
    
figure()
imshow(bundleImage)
title('imoverlay') ;
% 
% for i = 1:length(seg_list_new)
%     hold on;
% %     plot(seg_list_new{i}(1:end,2), seg_list_new{i}(1:end,1), 'r', 'linewidth', 2);
% %     q=3;
% %     p=1;  
% %     P = [seg_list_new{q}(1,2), seg_list_new{q}(1,1)];
% %     polyout = polybuffer(P,'points',r);
% %     plot(polyout)
% %     plot(11.0731, 794.5101, 'o')
% %     %plot([seg_list_new{p}(1,2) seg_list_new{p}(2,2)], [seg_list_new{p}(1,1) seg_list_new{p}(2,1)], 'g', 'linewidth', 2);
% %     plot([seg_list_new{q}(1,2) seg_list_new{q}(2,2)], [seg_list_new{q}(1,1) seg_list_new{q}(2,1)], 'g', 'linewidth', 2);
% %     plot([seg_list_new{p}(1,2) seg_list_new{p}(2,2)], [seg_list_new{p}(1,1) seg_list_new{p}(2,1)], 'g', 'linewidth', 2);
%     hold off;
% end

% %% To shapefile
% STR = 'struct(''Geometry'',values ,''X'', values,''Y'', values,''ID'',values)';
% values = cell(length(seg_list_new), 1);
% S = eval(STR);
% 
% for i=1:length(seg_list_new)
%     S(i).Geometry='Line';
%     S(i).ID=i;
%     [x,y]=pix2map(R2,seg_list_new{i}(1:end,1),seg_list_new{i}(1:end,2));
%     S(i).X=[x;NaN]';
%     S(i).Y=[y;NaN]';
% end
% shapewrite(S,'Road_network_2.shp');