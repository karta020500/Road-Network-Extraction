close all;
clear all;
%bundleImage=imread(imagepath);
[bundleImage,R] = geotiffread('NDVI_binary.tif');
R2 = worldfileread('NDVI_binary.tfw', 'planar', size(bundleImage));

% bundleImage(bundleImage > 3) = 1;
bundleImage(bundleImage == 65536) = 0;
[m,n] = size(bundleImage);

myfigure = figure('Name','Binary_ndvi_reflec.tif','NumberTitle', 'off');
bundleForPeakFinding = bundleImage;

%subplot(2,2,1); imshow(bundleImage);

bundleImage=double(bundleImage);
subplot(2,2,1); mapshow(bundleImage,R2);
imagesc(bundleImage);
axis image off

%pre-processing------------------------------------------------------------

% opening
se = strel('square',2);
opened = imopen(bundleImage,se);

%closeing
se2 = strel('square',5);
closed = imclose(opened,se2);

%connectivity_area
bundleImage = bwareaopen(closed, 5);
Pixel_connectivity = bundleImage - bwareaopen(closed, 3200);

bundleImage = Pixel_connectivity;

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
threshold_val = ceil(0.99*max(R(:)));

% initialize the loop variables
done = false;
h = R;
hnew = h;
nhood = size(R)/10;
nhood = max(2*ceil(nhood/2) + 1, 1);
nhood_center = (nhood-1)/2;
P = [];
numpeaks = 1;

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

peaks = [rho(P(:,1)) theta(P(:,2))];
plot(peaks(:,2), peaks(:,1), 's', 'color','white');
title('Radon Transform & Peaks') ;

subplot(2,2,4); plot(rho,R(:,P(1,2)));

[R_w,rho_w] = radon(bundleImage,peaks(2));


[p_m, p_n] = size(P);
[r_m, r_n] = size(R);
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

lines = radonlines(bundleImage,theta,rho,P,40,60,0);

for k = 1:numel(lines)
  seg_list{k}(1,1) = lines(k).point1(2); seg_list{k}(1,2) = lines(k).point1(1);
  seg_list{k}(2,1) = lines(k).point2(2); seg_list{k}(2,2) = lines(k).point2(1);
  %line([p1(1) p2(1)], [p1(2) p2(2)], 'Color', 'red', 'LineWidth', 2);
end
    
    
figure()
imshow(bundleImage)
title('imoverlay') ;
% 
for i = 1:length(seg_list)
    hold on;
    plot(seg_list{i}(1:end,2), seg_list{i}(1:end,1), 'r', 'linewidth', 2);
    hold off;
end

STR = 'struct(''Geometry'',values ,''X'', values,''Y'', values,''ID'',values)';
values = cell(length(seg_list), 1);
S = eval(STR);

for i=1:length(seg_list)
    S(i).Geometry='Line';
    S(i).ID=i;
    [x,y]=pix2map(R2,seg_list{i}(1:end,1),seg_list{i}(1:end,2));
    S(i).X=[x;NaN]';
    S(i).Y=[y;NaN]';
end
shapewrite(S,'try.shp');