function [closed, Width_seg]= GlidingBoxRadon(in)

[m,n] = size(in);

bundleForPeakFinding=in;
in=double(in);

% theta = [0:180]';
theta_step = 0.5;
theta = [0:theta_step:180-theta_step]';
[R, rho] = radon(in, theta);
%R = medfilt2(R,[3,3]);
R = imfilter(R, fspecial('gaussian', 3, 2));
%-------------------------------------------------------------

% Detect the peaks of transform output  
%Threshold value for peak detection
threshold_val = ceil(0.90*max(R(:))) ;
% Maximum nof peaks to identify in the image
done = false;
h = R;
hnew = h;
nhood = size(R)/5;
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

if  ~all(P == 0)

% Fast_P = FastPeakFind(R)
% 
% Px_index = Fast_P(1:2:end)
% Py_index = Fast_P(2:2:end)
% 
% peaks = [rho(Py_index) theta(Px_index)];

% Peak profile

    [p_m, p_n] = size(P);
    [r_m, r_n] = size(R);
    range = 5;
    semi_peak_radon = 0.9;
    Width_seg={};
    count=1;
    p_width=0;
% 
%linear interpolation
    for i = 1:p_m
        query = R(P(i,1),P(i,2))*semi_peak_radon;
        if P(i,1) > range && P(i,1) < r_m-range && P(i,2) > range && P(i,2) < r_n-range
        %if P(i,2) > range && P(i,2) < r_n-range

           p_range = R(P(i,1)-range:P(i,1)+range,P(i,2));
%            the_range = R(P(i,1),P(i,2)-range:P(i,2)+range);

           [~,p_range_I]=sort(abs(p_range-query));
            p_width = p_range(p_range_I(1)) - p_range(p_range_I(2));
%            [~,the_range_I]=sort(abs(the_range-query));

%            p_query = (p_range(p_range_I(1)) + p_range(p_range_I(2)))/2;
%            [~,p_query_I] = min(abs(p_range - p_query));
%            the_query = (the_range(the_range_I(1)) + the_range(the_range_I(2)))/2;
%            [~,the_query_I] = min(abs(the_range - the_query));

%            p0 = p_range(p_query_I);
%            p0_I = find(p_range == p0,1);
% 
%            the0 = the_range(the_query_I);
%            the0_I = find(the_range == the0,1);

           %P(i,1) = P(i,1) + ((range+1) - p0_I);
           %P(i,2) = P(i,2) + ((range+1) - the0_I);
           
           Width_seg{count}(1,1) = i;
           Width_seg{count}(1,2) = abs(p_width);
           count = count+1;
        end
    end

    peaks = [rho(P(:,1)) theta(P(:,2))];

    [R,rho] = radon(in,peaks(2));
    width1=fwhm(rho,R);

    stringout=sprintf('Width:%0.2f',abs(p_width));
    disp(stringout)

    f= width1;
    
%     if abs(p_width) > 7
%        figure()
%        imshow(in);
%        title('original') 
%     end

    %-------------------------------------------------------------

    % Detected lines on the image

    x_center = floor(size(in, 2)/2) ;
    y_center = floor(size(in, 1)/2) ;

    [p_m, p_n] = size(peaks);
    
    x1 = [];
    y1 = [];
    x2 = [];
    y2 = [];

    for p=1:p_m

        x_1 = [-x_center-100, x_center+100] ;
        y_1 = [0, 0] ;

        % Shift at first
        x_1_shifted = x_1 ;
        y_1_shifted = [y_1(1)-peaks(p,1), y_1(2)-peaks(p,1)] ;

        % Rotate 
        peaks(p,2) = 90 - peaks(p,2) ;
        t=peaks(p,2)*pi/180;
        rotation_mat = [ cos(t) -sin(t) ; sin(t) cos(t) ] ;
        x_y_rotated = rotation_mat*[x_1_shifted; y_1_shifted] ;
        x_rotated = x_y_rotated(1,:) ;
        y_rotated = x_y_rotated(2,:) ;
        
        x1(p) = x_rotated(1)+x_center;
        y1(p) = y_rotated(1)+y_center;
        x2(p) = x_rotated(2)+x_center;
        y2(p) = y_rotated(2)+y_center;
    end
    %-------------------------------------------------------------

    % Create a line from point 1 to point 2
    line_m = zeros(m,n);
    for l=1:p_m
        spacing = 1;
        numSamples = ceil(sqrt((x2(l)-x1(l))^2+(y2(l)-y1(l))^2) / spacing);
        x = linspace(x1(l), x2(l), numSamples);
        y = linspace(y1(l), y2(l), numSamples);
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
            line_m(finaly(p),finalx(p)) = in(finaly(p),finalx(p));
            end
        end
    end
        %se = strel('disk',9);
        %closed = imclose(line_m,se);
        closed = line_m;
        %Endpoint = bwmorph(closed,'Endpoints');
        %T = find(Endpoint==1);
        %final_L = [];
        %for i=1:length(T)
        %   F = fix(T(i)/m);
        %   L = abs(T(i)-(m*F));
        %   final_L(i,1) = L;
        %   final_L(i,2) = F;
        %end
        
    else
    closed = zeros(m,n);
end



