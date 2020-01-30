close all;
clear all;
[bundleImage,R] = geotiffread('Road_calss.tif');
R2 = worldfileread('Road_calss.tfw', 'planar', size(bundleImage));

bundleImage(bundleImage > 3) = 0;
bundleImage(bundleImage == 3) = 1;
bundleImage=double(bundleImage);

[ih, iw] = size(bundleImage);

%% pre-processing

%opening
se = strel('square',2);
opened = imopen(bundleImage,se);

%closeing
se2 = strel('square',5);
closed = imclose(opened,se2);
%closed = imclose(closed,se2);
%connectivity_area
Pixel_connectivity = bwareaopen(closed, 50);


%hole_fill
filled = imfill(Pixel_connectivity,'holes');
holes = filled & ~Pixel_connectivity;
bigholes = bwareaopen(holes, 200);
smallholes = holes & ~bigholes;
filled = Pixel_connectivity | smallholes;

skel = bwmorph(filled,'skel',Inf);

skel = bwmorph(skel, 'spur', 9);

skel = imdilate(skel,ones(5,5));

% figure()
% imshow(bundleImage);

%% GlidingBoxRadon(bundleImage);
m = 40;
n =  30;
block_image = skel;
fun = @(block_struct) GlidingBoxRadon(block_struct.data);
%block_radon = blockproc(block_image ,[m,n] ,fun);
block_radon = blockproc(block_image ,[m,n] ,fun, 'BorderSize', [10,10]);

se3 = strel('square',1);
block_radon = imclose(block_radon,se3);
block_radon = bwareaopen(block_radon, 20);

%% edgelinking
 block_radon = filledgegaps(block_radon, 3);
[edgelist, labelededgeim] = edgelink(block_radon, 10);

% Display the edgelists with random colours for each distinct edge 
% in figure 2

% drawedgelist(edgelist, size(block_radon), 1, 'rand', 2); axis off 

    % Fit line segments to the edgelists
tol = 10;         % Line segments are fitted with maximum deviation from
         % original edge of 2 pixels.
seglist = lineseg(edgelist, tol);

% Draw the fitted line segments stored in seglist in figure window 3 with
% a linewidth of 2 and random colours
% drawedgelist(seglist, size(block_radon), 2, 'rand', 3); axis off

%% disolve lines
polylines ={};
number = 1;
add_seg = {};
count = 1;

 for i=1:length(seglist)
     if ~isequal(size(seglist{i}),[2,2])
        polylines{number} = i;
        number = number+1;
     end
 end
 
 for i=1:length(polylines)
     seg = length(seglist{polylines{i}});
     for j=1:seg-2
         add_seg{count} = [seglist{polylines{i}}(j+2,:);seglist{polylines{i}}(j+1,:)];
         count=count+1;
     end
 end
 
for i = 1:length(add_seg)  
    seglist{length(seglist)+1} = add_seg{i};
end

%% Unique lines
% Strategy
currentfilename = 'Road_calss.tif';
[w, long_seg] = linewidthofimages(currentfilename);

seg_index = [];
for i = 1:length(long_seg)
    hold on;
    %plot([long_seg{i}(1,2) long_seg{i}(2,2)], [long_seg{i}(1,1) long_seg{i}(2,1)], 'r', 'linewidth', 1);
    for j = 1:length(seglist)
        midp = [(seglist{j}(1,1)+seglist{j}(2,1))/2,(seglist{j}(1,2)+seglist{j}(2,2))/2];
        x = midp; %some point
        a = long_seg{i}(1,:); %segment points a,b
        b = long_seg{i}(2,:);
        d_ab = norm(a-b);
        d_ax = norm(a-x);
        d_bx = norm(b-x);

        if dot(a-b,x-b)*dot(b-a,x-a)>=0
            A = [a,1;b,1;x,1];
            distance = abs(det(A))/d_ab;        
        else
            distance = min(d_ax, d_bx);
        end
%           distance = point_to_line_distance([midp(1,1),midp(1,2),0], [long_seg{i}(1,1),long_seg{i}(1,2),0],[long_seg{i}(2,1),long_seg{i}(2,2),0]);
            if distance <= 15
            diff_ang = (atan((long_seg{i}(2,1)-long_seg{i}(1,1))/(long_seg{i}(2,2)-long_seg{i}(1,2))) - atan((seglist{j}(2,1)-seglist{j}(1,1))/(seglist{j}(2,2)-seglist{j}(1,2)))) * 180/pi;
            if abs(diff_ang) >= 150 || abs(diff_ang) <= 30
               %removeing_index
                seg_index = [seg_index;j];
            end
            end
        %plot([seglist{j}(1,2) seglist{j}(2,2)], [seglist{j}(1,1) seglist{j}(2,1)], 'g', 'linewidth', 1);
    end
       l=9;
%      s=102;
%      plot([long_seg{l}(1,2) long_seg{l}(2,2)], [long_seg{l}(1,1) long_seg{l}(2,1)], 'g', 'linewidth', 1);
%      plot([seglist{s}(1,2) seglist{s}(2,2)], [seglist{s}(1,1) seglist{s}(2,1)], 'g', 'linewidth', 1);
     hold off;
end

seglist(seg_index) =[];

%% linkage
add_seg = {};
count = 1;

add_seg{count} = [366,135;370,60];%short
count =count+1;
% add_seg{count} = [805,770;816,680];%short
% count =count+1;
add_seg{count} = [840,370;820,366];%short
count =count+1;
% add_seg{count} = [126,128;29,133];%short
% count =count+1;
add_seg{count} = [450,1050;525,1050];%short
count =count+1;
add_seg{count} = [979,1325;980,1350];%short
count =count+1;
add_seg{count} = [985,100;923,92];%short
count =count+1;
add_seg{count} = [420,916;420,894];%short
count =count+1;
add_seg{count} = [854,411;845,480];%short
count =count+1;


for i = 1:length(add_seg)  
    seglist{length(seglist)+1} = add_seg{i};
 end

remove = [1123,731;1121,764];%short
index = cellfun(@(xx)isequal(xx,remove),seglist);
[row,col] = find(index);
if ~isempty(col)
   for i = 1:length(col)
     seglist(col(i)) = [];
   end
end

remove = [1123,811;1131,895];%short
index = cellfun(@(xx)isequal(xx,remove),seglist);
[row,col] = find(index);
if ~isempty(col)
   for i = 1:length(col)
     seglist(col(i)) = [];
   end
end

remove = [1123,781;1125,806];%short
index = cellfun(@(xx)isequal(xx,remove),seglist);
[row,col] = find(index);
if ~isempty(col)
   for i = 1:length(col)
     seglist(col(i)) = [];
   end
end

remove = [1127,914;1130,960];%short
index = cellfun(@(xx)isequal(xx,remove),seglist);
[row,col] = find(index);
if ~isempty(col)
   for i = 1:length(col)
     seglist(col(i)) = [];
   end
end

remove = [1138,1471;1127,1489];%short
index = cellfun(@(xx)isequal(xx,remove),seglist);
[row,col] = find(index);
if ~isempty(col)
   for i = 1:length(col)
     seglist(col(i)) = [];
   end
end

remove = [923,92;986,110];%short
index = cellfun(@(xx)isequal(xx,remove),seglist);
[row,col] = find(index);
if ~isempty(col)
   for i = 1:length(col)
     seglist(col(i)) = [];
   end
end

remove = [952,661;957,681];%short
index = cellfun(@(xx)isequal(xx,remove),seglist);
[row,col] = find(index);
if ~isempty(col)
   for i = 1:length(col)
     seglist(col(i)) = [];
   end
end

remove = [859,211;869,228];%short
index = cellfun(@(xx)isequal(xx,remove),seglist);
[row,col] = find(index);
if ~isempty(col)
   for i = 1:length(col)
     seglist(col(i)) = [];
   end
end

% remove = [861,391;854,411];%short
% index = cellfun(@(xx)isequal(xx,remove),seglist);
% [row,col] = find(index);
% if ~isempty(col)
%    for i = 1:length(col)
%      seglist(col(i)) = [];
%    end
% end

remove = [848,421;845,480];%short
index = cellfun(@(xx)isequal(xx,remove),seglist);
[row,col] = find(index);
if ~isempty(col)
   for i = 1:length(col)
     seglist(col(i)) = [];
   end
end

for i = 1:length(seglist)
    for j = 1:length(seglist)
         Distance = DistBetween2Segment(seglist{i}(1, :), seglist{i}(2, :), seglist{j}(1, :), seglist{j}(2, :));
%Distance between middle point of two line segments
         if abs(Distance)> 0 && abs(Distance) < 80
            diff_ang_seg = (atan((seglist{i}(2,1)-seglist{i}(1,1))/(seglist{i}(2,2)-seglist{i}(1,2))) - atan((seglist{j}(2,1)-seglist{j}(1,1))/(seglist{j}(2,2)-seglist{j}(1,2)))) * 180/pi;
            if abs(diff_ang_seg) >= 160 || abs(diff_ang_seg) <=23
               dist=bsxfun(@hypot,seglist{j}(:,1)-seglist{i}(1,1),seglist{j}(:,2)-seglist{i}(1,2));
               out = seglist{j}(dist==min(dist),:);
               dist=bsxfun(@hypot,seglist{i}(:,1)-out(1,1),seglist{i}(:,2)-out(1,2));
               final = seglist{i}(find(dist==min(dist), 1),:);
               diff_ang_conn_li = (atan((seglist{i}(2,1)-seglist{i}(1,1))/(seglist{i}(2,2)-seglist{i}(1,2))) - atan((out(1,1)-final(1,1))/(out(1,2)-final(1,2)))) * 180/pi;
               diff_ang_conn_lj = (atan((seglist{j}(2,1)-seglist{j}(1,1))/(seglist{j}(2,2)-seglist{j}(1,2))) - atan((out(1,1)-final(1,1))/(out(1,2)-final(1,2)))) * 180/pi;
%Lateral offset and direction difference of line segment
               if abs(Distance)> 0 && abs(Distance) < 5
                  if abs(d)> 0 && abs(d) < 150
                     add_seg{count}(1,:) = out;
                     add_seg{count}(2,:) = final;
                     count = count+1;
                  end
               elseif (abs(diff_ang_conn_li) >= 160 || abs(diff_ang_conn_li) <=20) && (abs(diff_ang_conn_lj) >= 160 || abs(diff_ang_conn_lj) <=20)
                  d = pdist([out(1,:); final(1,:)],'euclidean');
                  if abs(d)> 0 && abs(d) < 150
                     add_seg{count}(1,:) = out;
                     add_seg{count}(2,:) = final;
                     count = count+1;
                  end
               end
            end
        end
    end
end

add_seg{count} = [396,1440;401,1437];%main
count = count+1;

[add_seg,idx,idx2] = uniquecell(add_seg);
for i = 1:length(add_seg)  
    seglist{length(seglist)+1} = add_seg{i};
end

%% Merge line segments
% Create a line from point 1 to point 2
    line_m = zeros(ih,iw);
    
    % seglist
    for i=1:length(seglist)
        spacing = 0.5;
        numSamples = ceil(sqrt((seglist{i}(2,2)-seglist{i}(1,2))^2+(seglist{i}(2,1)-seglist{i}(1,1))^2) / spacing);
        x = linspace(seglist{i}(1,2), seglist{i}(2,2), numSamples);
        y = linspace(seglist{i}(1,1), seglist{i}(2,1), numSamples);
        xy = round([x',y']);
        dxy = abs(diff(xy, 1));
        duplicateRows = [0; sum(dxy, 2) == 0];
        % Now for the answer:
        finalxy = xy(~duplicateRows,:);
        finalx = finalxy(:, 1);
        finaly = finalxy(:, 2);
        % Plot the points

        for p=1:length(finalxy)
            if  finalxy(p,1) > 0 && finalxy(p,1) < iw && finalxy(p,2) > 0 && finalxy(p,2) < ih
            line_m(finaly(p),finalx(p)) = 1;
            end
        end
    end
    
    % long_seg
    for i=1:length(long_seg)
        spacing = 0.5;
        numSamples = ceil(sqrt((long_seg{i}(2,2)-long_seg{i}(1,2))^2+(long_seg{i}(2,1)-long_seg{i}(1,1))^2) / spacing);
        x = linspace(long_seg{i}(1,2), long_seg{i}(2,2), numSamples);
        y = linspace(long_seg{i}(1,1), long_seg{i}(2,1), numSamples);
        xy = round([x',y']);
        dxy = abs(diff(xy, 1));
        duplicateRows = [0; sum(dxy, 2) == 0];
        % Now for the answer:
        finalxy = xy(~duplicateRows,:);
        finalx = finalxy(:, 1);
        finaly = finalxy(:, 2);
        % Plot the points

        for p=1:length(finalxy)
            if  finalxy(p,1) > 0 && finalxy(p,1) < iw && finalxy(p,2) > 0 && finalxy(p,2) < ih
            line_m(finaly(p),finalx(p)) = 1;
            end
        end
    end
    
line_m = filledgegaps(line_m, 3);
[edgelist, labelededgeim] = edgelink(line_m, 1);

seg_list_new = lineseg(edgelist, 5);

%% linkage two
add_seg = {};
count = 1;
% 
% % for i = 1:length(seg_list_new)
% %     for j = 1:length(seg_list_new)
% %          Distance = DistBetween2Segment(seg_list_new{i}(1, :), seg_list_new{i}(2, :), seg_list_new{j}(1, :), seg_list_new{j}(2, :));
% % %Distance between middle point of two line segments
% %          if abs(Distance)> 0 && abs(Distance) < 80
% %             diff_ang_seg = (atan((seg_list_new{i}(2,1)-seg_list_new{i}(1,1))/(seg_list_new{i}(2,2)-seg_list_new{i}(1,2))) - atan((seg_list_new{j}(2,1)-seg_list_new{j}(1,1))/(seg_list_new{j}(2,2)-seg_list_new{j}(1,2)))) * 180/pi;
% %             if abs(diff_ang_seg) >= 160 || abs(diff_ang_seg) <=20
% %                dist=bsxfun(@hypot,seg_list_new{j}(:,1)-seg_list_new{i}(1,1),seg_list_new{j}(:,2)-seg_list_new{i}(1,2));
% %                out = seg_list_new{j}(dist==min(dist),:);
% %                dist=bsxfun(@hypot,seg_list_new{i}(:,1)-out(1,1),seg_list_new{i}(:,2)-out(1,2));
% %                final = seg_list_new{i}(find(dist==min(dist), 1),:);
% %                diff_ang_conn_lj = (atan((seg_list_new{j}(2,1)-seg_list_new{j}(1,1))/(seg_list_new{j}(2,2)-seg_list_new{j}(1,2))) - atan((out(1,1)-final(1,1))/(out(1,2)-final(1,2)))) * 180/pi;
% %                diff_ang_conn_li = (atan((seg_list_new{i}(2,1)-seg_list_new{i}(1,1))/(seg_list_new{i}(2,2)-seg_list_new{i}(1,2))) - atan((out(1,1)-final(1,1))/(out(1,2)-final(1,2)))) * 180/pi;
% % %Lateral offset and direction difference of line segment
% %                if (abs(diff_ang_conn_li) >= 170 || abs(diff_ang_conn_li) <=10) && (abs(diff_ang_conn_lj) >= 170 || abs(diff_ang_conn_lj) <=10)
% %                   d = pdist([out(1,:); final(1,:)],'euclidean');
% %                   if abs(d)> 0 && abs(d) < 150
% %                      add_seg{count}(1,:) = out;
% %                      add_seg{count}(2,:) = final;
% %                      count = count+1;
% %                   end
% %                end
% %             end
% %         end
% %     end
% % end
% 
add_seg{count} = [200,746;381,705];%main
count = count+1;
add_seg{count} = [569,741;479,718];%main
count = count+1;
add_seg{count} = [834,577;845,480];%main
count = count+1;
add_seg{count} = [861,301;863,210];%main
count = count+1;
% add_seg{count} = [676,508;706,386];%main
% count = count+1;
add_seg{count} = [713,361;750,203];%main
count = count+1;
% add_seg{count} = [797,270;860,284];%main
% count = count+1;
add_seg{count} = [980,1350;984,1415];%main
count = count+1;
add_seg{count} = [980,1212;979,1325];%main
count = count+1;
add_seg{count} = [1025,298;1030,256];%main
count = count+1;
add_seg{count} = [676,599;706,386];%main
count = count+1;
add_seg{count} = [801,789;807,744];%main
count = count+1;
add_seg{count} = [776,278;765,271];%main
count = count+1;


[add_seg,idx,idx2] = uniquecell(add_seg);
for i = 1:length(add_seg)  
    seg_list_new{length(seg_list_new)+1} = add_seg{i};
end

%% disolve lines
polylines ={};
number = 1;
add_seg = {};
count = 1;

 for i=1:length(seg_list_new)
     if ~isequal(size(seg_list_new{i}),[2,2])
        polylines{number} = i;
        number = number+1;
     end
 end
 
 for i=1:length(polylines)
     seg = length(seg_list_new{polylines{i}});
     for j=1:seg-2
         add_seg{count} = [seg_list_new{polylines{i}}(j+2,:);seg_list_new{polylines{i}}(j+1,:)];
         count=count+1;
     end
 end
 
for i = 1:length(add_seg)  
    seg_list_new{length(seg_list_new)+1} = add_seg{i};
end
 
%% Merge list all
% Create a line from point 1 to point 2
    line_m = zeros(ih,iw);
    
    % seglist
    for i=1:length(seg_list_new)
        spacing = 0.5;
        numSamples = ceil(sqrt((seg_list_new{i}(2,2)-seg_list_new{i}(1,2))^2+(seg_list_new{i}(2,1)-seg_list_new{i}(1,1))^2) / spacing);
        x = linspace(seg_list_new{i}(1,2), seg_list_new{i}(2,2), numSamples);
        y = linspace(seg_list_new{i}(1,1), seg_list_new{i}(2,1), numSamples);
        xy = round([x',y']);
        dxy = abs(diff(xy, 1));
        duplicateRows = [0; sum(dxy, 2) == 0];
        % Now for the answer:
        finalxy = xy(~duplicateRows,:);
        finalx = finalxy(:, 1);
        finaly = finalxy(:, 2);
        % Plot the points

        for p=1:length(finalxy)
            if  finalxy(p,1) > 0 && finalxy(p,1) < iw && finalxy(p,2) > 0 && finalxy(p,2) < ih
            line_m(finaly(p),finalx(p)) = 1;
            end
        end
    end
    
    line_m = filledgegaps(line_m, 3);
    [edgelist, labelededgeim] = edgelink(line_m, 1);

    seg_list_new = lineseg(edgelist, 5);
    
remove = [881,369;893,365];%short
index = cellfun(@(xx)isequal(xx,remove),seg_list_new);
[row,col] = find(index);
if ~isempty(col)
   for i = 1:length(col)
      seg_list_new(col(i)) = [];
   end
end
remove = [861,391;857,399];%short
index = cellfun(@(xx)isequal(xx,remove),seg_list_new);
[row,col] = find(index);
if ~isempty(col)
   for i = 1:length(col)
      seg_list_new(col(i)) = [];
   end
end

 remove = [861,391;857,396];%short
index = cellfun(@(xx)isequal(xx,remove),seg_list_new);
[row,col] = find(index);
if ~isempty(col)
   for i = 1:length(col)
      seg_list_new(col(i)) = [];
   end
end

%% remove short lines
% move_seg = {};
% count = 1;
%     for i=1:length(seg_list_new)
%         d = pdist([seg_list_new{i}(1,2), seg_list_new{i}(1,1);seg_list_new{i}(2,2), seg_list_new{i}(2,1)],'euclidean');
%         if abs(d)> 0 && abs(d) < 30
%            move_seg{count} = seg_list_new(i);
%            count = count +1;
%         end
%     end
%     
%     for i = 1:length(move_seg)  
%        remove = move_seg{i};
%        index = cellfun(@(xx)isequal(xx,remove{1,1}),seg_list_new);
%        [row,col] = find(index);
%        for j = 1:length(col)
%          seg_list_new(col(j)) = [];
%        end
%     end

%% Junction connection
% r = 55;
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
%    %                  [checkPt_i, onEnd] = checkPointOnSegment(segment_i,point,0);
%    %                  [checkPt_j, onEnd] = checkPointOnSegment(segment_j,point,0);
%                     polyout_2 = polybuffer(point,'points',1);
%                     [checkPt_i,out] = intersect(polyout_2,segment_i);
%                     [checkPt_j,out] = intersect(polyout_2,segment_j);
%                     %angle = (atan((seg_list_new{j}(2,1)-seg_list_new{j}(1,1))/(seg_list_new{j}(2,2)-seg_list_new{j}(1,2))) - atan((seg_list_new{j}(2,1)-seg_list_new{j}(1,1))/(seg_list_new{j}(2,2)-seg_list_new{j}(1,2)))) * 180/pi;
%                     disten = pdist([P(1,:); point(1,:)],'euclidean');
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
% %                  elseif abs(disten) <= 20
% %                        dist=bsxfun(@hypot,seg_list_new{i}(:,2)-point(1,1),seg_list_new{i}(:,1)-point(1,2));
% %                         if dist(1) == dist(2)
% %                            dist = unique(dist);
% %                         end
% %                        short_end = seg_list_new{i}(dist==min(dist),:);
% %                        %seg_list_new{i}(dist==min(dist),:) = point;
% %                        d = pdist([short_end(1,2), short_end(1,1); point(1,:)],'euclidean');
% %                        if abs(d)> 0 && abs(d) < r
% %                           add_seg{count}(1,:) = short_end;
% %                           add_seg{count}(2,:) = [point(1,2),point(1,1)];
% %                           count = count+1;
% %                        end
% %                        dist=bsxfun(@hypot,seg_list_new{j}(:,2)-point(1,1),seg_list_new{j}(:,1)-point(1,2));
% %                         if dist(1) == dist(2)
% %                            dist = unique(dist);
% %                         end
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
%              P_2 = [seg_list_new{i}(end,2), seg_list_new{i}(end,1)];
%              line_seg = [seg_list_new{j}(1:end,2) seg_list_new{j}(1:end,1)];
%              polyout = polybuffer(P_2,'points',r);
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
%                     %angle = (atan((seg_list_new{j}(2,1)-seg_list_new{j}(1,1))/(seg_list_new{j}(2,2)-seg_list_new{j}(1,2))) - atan((seg_list_new{j}(2,1)-seg_list_new{j}(1,1))/(seg_list_new{j}(2,2)-seg_list_new{j}(1,2)))) * 180/pi;
%                     disten = pdist([P_2(1,:); point(1,:)],'euclidean');
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
% %                  elseif abs(disten) <= 20
% %                        dist=bsxfun(@hypot,seg_list_new{i}(:,2)-point(1,1),seg_list_new{i}(:,1)-point(1,2));
% %                         if dist(1) == dist(2)
% %                            dist = unique(dist);
% %                         end
% %                        short_end = seg_list_new{i}(dist==min(dist),:);
% %                        %seg_list_new{i}(dist==min(dist),:) = point;
% %                        d = pdist([short_end(1,2), short_end(1,1); point(1,:)],'euclidean');
% %                        if abs(d)> 0 && abs(d) < r
% %                           add_seg{count}(1,:) = short_end;
% %                           add_seg{count}(2,:) = [point(1,2),point(1,1)];
% %                           count = count+1;
% %                        end
% %                        dist=bsxfun(@hypot,seg_list_new{j}(:,2)-point(1,1),seg_list_new{j}(:,1)-point(1,2));
% %                         if dist(1) == dist(2)
% %                            dist = unique(dist);
% %                         end
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
% 
% add_seg{count} = [963,422;956,396];%else
% count = count+1;
% add_seg{count} = [963,422;975,420];%else
% count = count+1;
% add_seg{count} = [909,1017;904,1019];%else
% count = count+1;
% add_seg{count} = [909,1017;910,1051];%else
% count = count+1;
% add_seg{count} = [302,82;306,45];
% count =count+1;
% add_seg{count} = [121,592;106,596];
% count =count+1;
% add_seg{count} = [868,61;868,32];
% count =count+1;
% add_seg{count} = [28,1165;11,1192];
% count =count+1;
% add_seg{count} = [217,32;217,47];
% count =count+1;
% add_seg{count} = [217,32;200,32];
% count =count+1;
% add_seg{count} = [721,1346;676,1350];
% count =count+1;
% add_seg{count} = [676,1389;676,1456];
% count =count+1;
% add_seg{count} = [54,329;58,388];
% count =count+1;
% add_seg{count} = [326,1299;270,1299];
% count =count+1;
% add_seg{count} = [593,1321;593,1262];
% count =count+1;
% add_seg{count} = [972,734;1001,759];
% count =count+1;
% add_seg{count} = [257,1281;238,1267];
% count =count+1;
% add_seg{count} = [377,1355;326,1355];
% count =count+1;
% 
% [add_seg,idx,idx2] = uniquecell(add_seg);
% for i = 1:length(add_seg)  
%      seg_list_new{length(seg_list_new)+1} = add_seg{i};
% end
% 
% remove = [1074,499;1103,469];
% index = cellfun(@(xx)isequal(xx,remove),seg_list_new);
% [row,col] = find(index);
% if ~isempty(col)
%    for i = 1:length(col)
%       seg_list_new(col(i)) = [];
%    end
% end
% 
% remove = [881,369;950,355];
% index = cellfun(@(xx)isequal(xx,remove),seg_list_new);
% [row,col] = find(index);
% if ~isempty(col)
%    for i = 1:length(col)
%       seg_list_new(col(i)) = [];
%    end
% end
% 
% remove = [847,31;868,49];
% index = cellfun(@(xx)isequal(xx,remove),seg_list_new);
% [row,col] = find(index);
% if ~isempty(col)
%    for i = 1:length(col)
%       seg_list_new(col(i)) = [];
%    end
% end
% 
% remove = [450,1458;470,1455];
% index = cellfun(@(xx)isequal(xx,remove),seg_list_new);
% [row,col] = find(index);
% if ~isempty(col)
%    for i = 1:length(col)
%       seg_list_new(col(i)) = [];
%    end
% end
% 
% remove = [950,355;893,365];
% index = cellfun(@(xx)isequal(xx,remove),seg_list_new);
% [row,col] = find(index);
% if ~isempty(col)
%    for i = 1:length(col)
%       seg_list_new(col(i)) = [];
%    end
% end
% 
% remove = [881,369;893,365];
% index = cellfun(@(xx)isequal(xx,remove),seg_list_new);
% [row,col] = find(index);
% if ~isempty(col)
%    for i = 1:length(col)
%       seg_list_new(col(i)) = [];
%    end
% end
% 
% remove = [801,789;800,737];
% index = cellfun(@(xx)isequal(xx,remove),seg_list_new);
% [row,col] = find(index);
% if ~isempty(col)
%    for i = 1:length(col)
%       seg_list_new(col(i)) = [];
%    end
% end
% 
% remove = [1001,759;1023,732];
% index = cellfun(@(xx)isequal(xx,remove),seg_list_new);
% [row,col] = find(index);
% if ~isempty(col)
%    for i = 1:length(col)
%       seg_list_new(col(i)) = [];
%    end
% end
% remove = [1001,759;1023,732];
% index = cellfun(@(xx)isequal(xx,remove),seg_list_new);
% [row,col] = find(index);
% if ~isempty(col)
%    for i = 1:length(col)
%       seg_list_new(col(i)) = [];
%    end
% end
% 

%% adding
% add_seg={};
% count=1;
% 
% % add_seg{count} = [366,143;370,45];
% % count =count+1;
% % add_seg{count} = [801,789;807,744];
% % count =count+1;
% % add_seg{count} = [34,1226;61,1242];
% % count =count+1;
% % add_seg{count} = [794,270;861,284];
% % count =count+1;
% % add_seg{count} = [775,350;860,373];
% % count =count+1;
% % add_seg{count} = [1030,256;1025,298];
% % count =count+1;
% % add_seg{count} = [863,210;861,284];
% % count =count+1;
% % add_seg{count} = [980,1253;984,1415];
% % count =count+1;
% % add_seg{count} = [1080,1179;1134,1179];
% % count =count+1;
% % add_seg{count} = [843,504;834,577];
% % count =count+1;
% % add_seg{count} = [424,1050;546,1050];
% % count =count+1;
% % % add_seg{count} = [126,128;177,131];
% % % count =count+1;
% % add_seg{count} = [322,940;322,991];
% % count =count+1;
% % add_seg{count} = [419,938;419,989];
% % count =count+1;
% add_seg{count} = [54,329;58,388];
% count =count+1;
% % % add_seg{count} = [862,739;806,739];
% % % count =count+1;
% % % add_seg{count} = [920,739;972,739];
% % % count =count+1;
% add_seg{count} = [326,1299;270,1299];
% count =count+1;
% % add_seg{count} = [400,707;453,713];
% % count =count+1;
% % add_seg{count} = [478,718;538,741];
% % count =count+1;
% add_seg{count} = [593,1321;593,1262];
% count =count+1;
% add_seg{count} = [972,734;1001,759];
% count =count+1;
% add_seg{count} = [257,1281;238,1267];
% count =count+1;
% % add_seg{count} = [753,1343;805,1340];
% % count =count+1;
% add_seg{count} = [377,1355;326,1355];
% count =count+1;
% % add_seg{count} = [94,775;170,752];
% % count =count+1;
% % add_seg{count} = [228,740;381,705];
% % count =count+1;
% 
% 
% for i = 1:length(add_seg)  
%      seg_list_new{length(seg_list_new)+1} = add_seg{i};
% end

%% removeing
% remove = [859,211;835,171];
% index = cellfun(@(xx)isequal(xx,remove),seg_list_new);
% [row,col] = find(index);
% if ~isempty(col)
%    for i = 1:length(col)
%      seg_list_new{col(i)} = [0,0;0,0];
%    end
% end
% 
% remove = [859,211;869,228];
% index = cellfun(@(xx)isequal(xx,remove),seg_list_new);
% [row,col] = find(index);
% if ~isempty(col)
%    for i = 1:length(col)
%       seg_list_new{col(i)} = [0,0;0,0];
%    end
% end
% 
% remove = [1137,1473;1152,1446];
% index = cellfun(@(xx)isequal(xx,remove),seg_list_new);
% [row,col] = find(index);
% if ~isempty(col)
%    for i = 1:length(col)
%       seg_list_new{col(i)} = [0,0;0,0];
%    end
% end
% 
% remove = [1128,1488;1137,1473];
% index = cellfun(@(xx)isequal(xx,remove),seg_list_new);
% [row,col] = find(index);
% if ~isempty(col)
%    for i = 1:length(col)
%       seg_list_new{col(i)} = [0,0;0,0];
%    end
% end
% 
% remove = [1121,499;1074,500];
% index = cellfun(@(xx)isequal(xx,remove),seg_list_new);
% [row,col] = find(index);
% if ~isempty(col)
%    for i = 1:length(col)
%       seg_list_new{col(i)} = [0,0;0,0];
%    end
% end
% 
% remove = [126,128;177,131];
% index = cellfun(@(xx)isequal(xx,remove),seg_list_new);
% [row,col] = find(index);
% if ~isempty(col)
%    for i = 1:length(col)
%     seg_list_new{col(i)} = [0,0;0,0];
%    end
% end

%% drawedlines

% figure()
% imshow(filled)
% title('Result') ;

for i = 1:length(seg_list_new)
    hold on;
    plot(seg_list_new{i}(1:end,2), seg_list_new{i}(1:end,1), 'g', 'linewidth', 2);
%     plot([731, 741], [1153 ,538], 'r', 'linewidth', 2);
%   q=125;
%   k=200;
%   P = [seg_list_new{q}(2,2), seg_list_new{q}(2,1)];
%   polyout = polybuffer(P,'points',r);
%   plot(seg_list_new{q}(2,2), seg_list_new{q}(2,1), 'o')
%   plot(1345, 1152, 'o')
%   plot(polyout)
%   plot([seg_list_new{q}(1,2) seg_list_new{q}(2,2)], [seg_list_new{q}(1,1) seg_list_new{q}(2,1)], 'b', 'linewidth', 2);
%   plot([seg_list_new{k}(1,2) seg_list_new{k}(2,2)], [seg_list_new{k}(1,1) seg_list_new{k}(2,1)], 'r', 'linewidth', 2);
    hold off;
end

% for i = 1:length(long_seg)
%     hold on;
%     plot(long_seg{i}(1:end,2), long_seg{i}(1:end,1), 'r', 'linewidth', 3);
%     hold off;
% end

%%  figure-
% drawedgelist(seg_list_new, size(line_m), 2, 'rand', 3); axis off

figure()
subplot(2,2,1);
mapshow(bundleImage,R2);
imagesc(bundleImage);
title('original') ;

subplot(2,2,2);
mapshow(opened,R2);
imagesc(opened);
title('opening') ;

subplot(2,2,3);
mapshow(closed,R2);
imagesc(closed);
title('closeing') ;

subplot(2,2,4);
mapshow(filled,R2);
imagesc(filled);
title('hole filling') ;


figure()
subplot(1,2,1);
overlay = imoverlay(filled, block_radon, 'r');
mapshow(overlay, R2);
imagesc(overlay);
title('input image with center line') ;

% subplot(1,2,2);
% mapshow(line_m, R2);
% imagesc(line_m);
% title('center line') ;

% figure()
% mapshow(skel, R2);
% imagesc(skel);
% title('skel') ;

figure()
img = zeros(1171,1489);
img(m:m:end,:,:) = 255;       %# Change every tenth row to black
img(:,n:n:end,:) = 255;       %# Change every tenth column to black
overlay_2 = imoverlay(overlay, img, 'g');
mapshow(overlay_2, R2);
imagesc(overlay_2);
title('grid prosessing') ;
% 
% figure()
% overlay_2 = imoverlay(skel, img, 'g');
% mapshow(overlay_2, R2);
% imagesc(overlay_2);
% title('grid prosessing') ;

%% To shapefile
STR = 'struct(''Geometry'',values ,''X'', values,''Y'', values,''ID'',values)';
values = cell(length(seg_list_new), 1);
S = eval(STR);

for i=1:length(seg_list_new)
    S(i).Geometry='Line';
    S(i).ID=i;
    [x,y]=pix2map(R2,seg_list_new{i}(1:end,1),seg_list_new{i}(1:end,2));
    S(i).X=[x;NaN]';
    S(i).Y=[y;NaN]';
end
shapewrite(S,'Road_network_2.shp');
