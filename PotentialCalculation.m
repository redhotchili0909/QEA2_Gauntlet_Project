clear
clf

%Potential Contour Plot
x_coord = -1.200:0.1:1.207;
y_coord = -1.213:0.1:1.092;
j=1;

for  x = -1.200:0.1:1.207
    i=1;
    for y = -1.213:0.1:1.092
        [pot(i,j),~]=pointDetail(x,y);
        i = i +1;
    end
    j = j +1;
end
contour(x_coord,y_coord,pot,20)

%Gradient Quiver Plot
[gradx_test,grady_test] = gradient(pot);
quiver(x_coord,y_coord,gradx_test,grady_test)

%Gradient Decent Calculation
delta = 0.41;
lambda = 0.09;
tolerance = 0.01;
n_max = 20;

r_i = [-0.838;-0.933];
points = r_i;
n= 0;

[~,grad_next] = pointDetail(0,0);
r_i = r_i + (lambda .* grad_next);

while(n < n_max && (norm(grad_next) > tolerance))
    r_i = r_i + (lambda .* grad_next);
    points = [points,r_i];
    [~,grad_next] = pointDetail(r_i(1),r_i(2));
    lambda = delta * lambda;
    n = n + 1;
end

hold on
contour(x_coord,y_coord,pot,20)
plot(points(1,:),points(2,:),".","Color","red")

clf
clear

%Neato Navigation Code

%[sensors,vels]=neatoSim(-0.838,-0.933,pi);
%[sensros,vels]=neato()

pause(10) 
delta = 0.41;
lambda = 0.09;
tolerance = 0.01;
n_max = 20;
r_i = [-0.838;-0.933];
n= 0;
grad_origin = [-1; 0];

[~,grad_next] = pointDetail(0,0);

find_theta = @(u,v) [atan2(dot((cross([u(1) u(2) 0], [v(1) v(2) 0])), [0 0 1]), dot([u(1) u(2) 0], [v(1) v(2) 0]))];

angle = find_theta(grad_origin, grad_next);

r_i = r_i + (lambda .* grad_next);

while(n < n_max && (norm(grad_next) > tolerance))
    angle = find_theta(grad_origin, grad_next);
    
    if (angle >= 0) % if angle is positive, turn counter-clockwise.
        vels.lrWheelVelocitiesInMetersPerSecond=[-0.05, 0.05];
    else % if angle is negative, turn clockwise
        vels.lrWheelVelocitiesInMetersPerSecond=[0.05, -0.05];
    end
    pause(abs(angle)/0.4081632653)

    vels.lrWheelVelocitiesInMetersPerSecond=[0.25, 0.25];
    pause(norm(lambda .* grad_next)/0.3)
    
    vels.lrWheelVelocitiesInMetersPerSecond=[0,0];
    pause(1)

    r_i = r_i + (lambda .* grad_next);
    
    grad_origin = grad_next;
    [~,grad_next] = pointDetail(r_i(1),r_i(2));
    lambda = delta * lambda;
    n = n + 1;
end

function [m,b] = slope_intercept(p1,p2) % P1 and P2 are Two sets of X and Y
    m = (p1(1)-p2(1))./(p1(2)-p2(2));
    b = p1(2)-m*p1(1);
end

function [pot,grad] = sloped_linePotential(strength,slope,intercept,x,y,t1,t2)
    fun1 = @(t,x,y,strength) strength.*log(sqrt((x-(t)).^2+(y-(intercept + slope*t)).^2))*sqrt((1).^2+(slope).^2);
    fun2 = @(t,x,y,strength) strength.*((x-t)./((x-t).^2+(y-(intercept + slope*t)).^2))*sqrt((1).^2+(slope).^2);
    fun3 = @(t,x,y,strength) strength.*((y-(intercept + slope*t))./((x-t).^2+(y-(intercept + slope*t)).^2))*sqrt((1).^2+(slope).^2);
    pot = integral(@(t) fun1(t,x,y,strength),t1,t2);
    grad = [integral(@(t) fun2(t,x,y,strength),t1,t2);integral(@(t) fun3(t,x,y,strength),t1,t2)];
end

function [pot,grad] = vertical_linePotential(strength,intercept,x,y,t1,t2)
    fun1 = @(t,x,y,strength) strength.*log(sqrt((x-intercept).^2+(y-(t)).^2))*1;
    fun2 = @(t,x,y,strength) strength.*((x-intercept)./((x-intercept).^2+(y-(t)).^2))*1;
    fun3 = @(t,x,y,strength) strength.*((y-t)./((x-intercept).^2+(y-(t)).^2))*1;
    pot = integral(@(t) fun1(t,x,y,strength),t1,t2);
    grad = [integral(@(t) fun2(t,x,y,strength),t1,t2);integral(@(t) fun3(t,x,y,strength),t1,t2)];
end

function [pot,grad] = circlePotential(strength,radius,x,y,center_p1,center_p2)
    fun1 = @(t,x,y,strength,radius) strength.*log(sqrt((x-(radius*cos(t)+center_p1)).^2+(y-(radius*sin(t)+center_p2)).^2))*radius;
    fun2 = @(t,x,y,strength) strength.*((x-(radius*cos(t)+center_p1))./((x-(radius*cos(t)+center_p1)).^2+(y-(radius*sin(t)+center_p2)).^2))*radius;
    fun3 = @(t,x,y,strength) strength.*((y-(radius*sin(t)+center_p2))./((x-(radius*cos(t)+center_p1)).^2+(y-(radius*sin(t)+center_p2)).^2))*radius;
    pot = integral(@(t) fun1(t,x,y,strength,radius),0,2*pi);
    grad = [integral(@(t) fun2(t,x,y,strength),0,2*pi);integral(@(t) fun3(t,x,y,strength),0,2*pi)];
end

function [pot,gradient] = pointDetail(x,y)
        object_strength = 1;
        wall_strength = 2;
        circle_strength = 100;

        %Surronding Walls and BoB
        [pot1,grad1]=vertical_linePotential(wall_strength,-1.200,x,y,-1.213,1.092);
        [pot2,grad2]=vertical_linePotential(wall_strength,1.207,x,y,-1.213,1.092);
        [pot3,grad3]=sloped_linePotential(wall_strength,0,1.092,x,y,-1.200,1.207);
        [pot4,grad4]=sloped_linePotential(wall_strength,0,-1.213,x,y,-1.200,1.207);
        [pot5,grad5]=circlePotential(circle_strength,0.01,x,y,1.031,-0.562);
        
        % Top left object
        top_left_1 = [-0.673,0.190];
        top_left_2 = [-0.255,0.609];
        top_left_3 = [-0.768,0.322];
        top_left_4 = [-0.429,0.735];
            % Bottom line
        [bottom_line_m,bottom_line_b] = slope_intercept(top_left_1,top_left_2);
        [pot6,grad6]=sloped_linePotential(object_strength,bottom_line_m,bottom_line_b,x,y,-0.6732,-0.2553);
            % Top line
        [top_line_m,top_line_b] = slope_intercept(top_left_3,top_left_4);
        [pot7,grad7]=sloped_linePotential(object_strength,top_line_m,top_line_b,x,y,-0.768,-0.429);
            % Left line
        [left_line_m,left_line_b] = slope_intercept(top_left_1,top_left_3);
        [pot8,grad8]=sloped_linePotential(object_strength,left_line_m,left_line_b,x,y,-0.673,-0.768);
            % Right line
        [right_line_m,right_line_b] = slope_intercept(top_left_4,top_left_2);
        [pot9,grad9] = sloped_linePotential(object_strength,right_line_m,right_line_b,x,y,-0.429,-0.255);
      
        % Bottom left object
        bottom_left_1 = [-0.151,-0.580];
        bottom_left_2 = [-0.013,-0.416];
        bottom_left_3 = [-0.374,-0.352];
        bottom_left_4 = [-0.211,-0.248];        
             % Bottom line
        [bottom_line_m,bottom_line_b] = slope_intercept(bottom_left_1,bottom_left_2);        
        [pot10,grad10]=sloped_linePotential(object_strength,bottom_line_m,bottom_line_b,x,y,-0.151, -0.013);
            % Top line
        [top_line_m,top_line_b] = slope_intercept(bottom_left_3,bottom_left_4);                 
        [pot11,grad11]=sloped_linePotential(object_strength,top_line_m,top_line_b,x,y,-0.374,-0.211);
            % Left line
        [left_line_m,left_line_b] = slope_intercept(bottom_left_1,bottom_left_3);
        [pot12,grad12]=sloped_linePotential(object_strength,left_line_m,left_line_b,x,y,-0.374, -0.151);
            % Right line
        [right_line_m,right_line_b] = slope_intercept(bottom_left_4,bottom_left_2);
        [pot13,grad13] = sloped_linePotential(object_strength,right_line_m,right_line_b,x,y,-0.211,-0.013);


        % Bottom right object
        bottom_right_1 = [0.458,-0.758];
        bottom_right_2 = [0.687,-0.524];
        bottom_right_3 = [0.288,-0.493];
        bottom_right_4 = [0.568,-0.363];        
             % Bottom line
        [bottom_line_m,bottom_line_b] = slope_intercept(bottom_right_1,bottom_right_2);        
        [pot14,grad14]=sloped_linePotential(object_strength,bottom_line_m,bottom_line_b,x,y,0.458, 0.687);
            % Top line
        [top_line_m,top_line_b] = slope_intercept(bottom_right_3,bottom_right_4);                 
        [pot15,grad15]=sloped_linePotential(object_strength,top_line_m,top_line_b,x,y,0.288,0.568);
            % Left line
        [left_line_m,left_line_b] = slope_intercept(bottom_right_1,bottom_right_3);
        [pot16,grad16]=sloped_linePotential(object_strength,left_line_m,left_line_b,x,y,0.458, 0.288);
            % Right line
        [right_line_m,right_line_b] = slope_intercept(bottom_right_4,bottom_right_2);
        [pot17,grad17] = sloped_linePotential(object_strength,right_line_m,right_line_b,x,y,0.568,0.687);


        % Top right object
        top_right_1 = [0.329,0.220];
        top_right_2 = [0.498,0.370];
        top_right_3 = [0.157,0.384];
        top_right_4 = [0.328,0.547];        
             % Bottom line
        [bottom_line_m,bottom_line_b] = slope_intercept(top_right_1,top_right_2);        
        [pot18,grad18]=sloped_linePotential(object_strength,bottom_line_m,bottom_line_b,x,y,0.458, 0.687);
            % Top line
        [top_line_m,top_line_b] = slope_intercept(top_right_3,top_right_4);                 
        [pot19,grad19]=sloped_linePotential(object_strength,top_line_m,top_line_b,x,y,0.288,0.568);
            % Left line
        [left_line_m,left_line_b] = slope_intercept(top_right_1,top_right_3);
        [pot20,grad20]=sloped_linePotential(object_strength,left_line_m,left_line_b,x,y,0.458, 0.288);
            % Right line
        [right_line_m,right_line_b] = slope_intercept(top_right_4,top_right_2);
        [pot21,grad21] = sloped_linePotential(object_strength,right_line_m,right_line_b,x,y,0.568,0.687);

        % Total Potential/Gradient
        pot = pot1+pot2+pot3+pot4-pot5+pot6+pot7+pot8+pot9+pot10+pot11+pot12+pot13+pot14+pot15+pot16+pot17+pot18+pot19+pot20+pot21;
        gradient =grad1+grad2+grad3+grad4-grad5+grad6+grad7+grad8+grad9+grad10+grad11+grad12+grad13+grad14+grad15+grad16+grad17+grad18+grad19+grad20+grad21;
end

