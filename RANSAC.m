clear
clf
load("finalgaunlet_map.mat");

max_round = 300;
max_count = 300;
count_test = 0;
thresh = 0.01; 
gap_thresh = 0.6;
xplot = all_points(1,:);
yplot = all_points(2,:);
count = 0;

hold off
plot(xplot,yplot,".","Color","red")
[radians,range] = cart2pol(xplot,yplot);
data = [radians;range];
[x,y]=pol2cart(data(1,:),data(2,:));
total=[x;y];    
exclude_line = 1;
while exclude_line
    [coeff,endp,inliners,outliers,frame]=robustLineFit(data,thresh,max_round);
    exclude_line = filter_out_gaps(inliners,1);
end
inliners_all = inliners;
endp_all = [];
[theta_re,rho_re] = cart2pol(outliers(1,:),outliers(2,:));
 new_data = [theta_re;rho_re];
 
while  count<max_count && length(new_data(1,:)) >1
    [coeff,endp,inliners,new_outlier,frame]=robustLineFit(new_data,thresh,max_round);
    if ~filter_out_gaps(inliners,gap_thresh)
        [theta_re,rho_re] = cart2pol(new_outlier(1,:),new_outlier(2,:));
        new_data = [theta_re;rho_re];
        inliners_all = [inliners_all,inliners];
        endp_all = [endp_all,endp];
        x1 = endp(1,1);
        y1 = endp(2,1);
        x2 = endp(1,2);
        y2 = endp(2,2);
        line([x1,x2],[y1,y2],'Color','green')
    end
    count = count + 1;
end

plot(inliners_all(1,:),inliners_all(2,:),".","Color","blue")


function [coeff,endp,inliners,outliers,frame]= robustLineFit(data,thresh,max_round)
   max_count = 0;
    pdistance = [];
    data(:,data(2,:)==0) = [];
    [x,y]=pol2cart(data(1,:),data(2,:));
    total= [x;y];
    
    for i = 1:max_round
        [v1,v2,index]=generate_points(x,y);
        for j= 1:length(total)
            pdistance(j)=point_to_line([total(:,j);0],[v1;0],[v2;0]);
        end
        inliner_count = sum(pdistance <= thresh);
        if inliner_count >= max_count
            max_count = inliner_count;
            max_index = index;
            max_inliner = pdistance <= thresh;
        end
    end
    
    p1 = [x(max_index(1));y(max_index(1))];
    p2 = [x(max_index(2));y(max_index(2))];
    diff = p2-p1;
    m = diff(1)./diff(2);
    frame = [p1,p2];
    coeff = [m;p1(2)-m*p1(1)];
    inliners = total(:,max_inliner);
    outliers = total(:,~max_inliner);
    [~,minx_index] = min(inliners(1,:));
    [~,maxx_index] = max(inliners(1,:));
    endp = [inliners(:,minx_index),inliners(:,maxx_index)];
end

function exclude_line = filter_out_gaps(inliers_all, gap_threshold)
    [~, order] = sort(inliers_all(1,:));
    inliers_all = inliers_all(:,order);

    inliers_dist_between_x = diff(inliers_all(1, :));
    inliners_dist_between_y = diff(inliers_all(2,:));
    if any(inliers_dist_between_x > gap_threshold) || any(inliners_dist_between_y >gap_threshold)
        exclude_line = true;
    else
        exclude_line = false;
    end
end
function [v1,v2,index] = generate_points(x,y)
    index = randi(length(x),[1,2]);
    v1 = [x(index(1));y(index(1))];
    v2 = [x(index(2));y(index(2))];
end

function distance = point_to_line(pt, v1, v2)
      a = v1 - v2;
      b = pt - v2;
      distance = norm(cross(a,b)) / norm(a);
end