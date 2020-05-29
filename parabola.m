% Function to determine to determine ball's trajectory and required
% velocities and torques depending on its position in the 3D space
%
% Input parameters:
%   p1 = coordinates in 3D space representing ball's center of gravity
%   p2 = coordinates in 3D space representign basket's location
%   theta = desired angle between ball and basket (end of trajectory)
%   m = mass of the object to be thrown 
%
% Output parameters:
%    v_rot = array of rotated vectors.


function [v, w] = parabola(p1, p2, theta, m)
    
    %Checks if there is any error in the input parameters
    if(length(p1)~=3)
        error('input vector is not three dimensional'), 
    end
    if(length(p2)~=3)
        error('input vector is not three dimensional'), 
    end
    
    %Translates the points so that p1 is in origin
    p0 = p1';
    p2 = p2 - p1';
    p1 = p1' - p1';
    
    %Determines the angle the coordinate system has to rotate in y axis for
    %the xy plane to be aligned with the direction of the shoot
    u = [1, 0, 0];
    w = [p2(1), 0, p2(3)];
    alpha = atan2(norm(cross(w,u)),dot(w,u));
      
    %Rotates the direction of the shot so that it is in XY plane
    temp = p2';
    T = [cos(alpha) 0 sin(alpha); 0 1 0; -sin(alpha) 0 cos(alpha)];
    p2 = T*temp;    
    if(abs(p2(3))>=1e-6)
        alpha = -alpha;
        T = [cos(alpha) 0 sin(alpha); 0 1 0; -sin(alpha) 0 cos(alpha)];
        p2 = T*temp;
    end
    
    %Determines the parabola parameters: R(1)x^2 + R(2)x + R(3)
    A = [(p1(1)^2) p1(1) 1; (p2(1)^2) p2(1) 1; 2*p2(1) 1 0];
    b = [p1(2); p2(2); -tan(theta)];
    R = inv(A)*b;
    
    %Obtains the highest point in the parabola
    delta = (R(2)^2) - 4*R(1)*R(3);
    ymax = -delta/(4*R(1));
    xmax = -R(2)/(2*R(1));
    max = [xmax; ymax; 0];
    
    %delta = (p2(1) - p1(1))/100;
    %for i=0:100
    %    x(i+1)=delta*i;
    %    y(i+1)=R(1)*(x(i+1)^2) + R(2)*x(i+1) + R(3);
    %end
    
    %Uses the maximal point to obtain the necessary speed:
    %The kinetic energy in the shooting moment must be equal to the
    %potential energy at the highest point
    g = 9.81;
    beta = atan(R(2));
    v0 = (sqrt(2*g*ymax))/sin(beta);

    %Rotates and translates the speed vector
    T = [cos(-alpha) 0 sin(-alpha); 0 1 0; -sin(-alpha) 0 cos(-alpha)];
    v = v0.*(T*[cos(beta); sin(beta); 0]);
    w = [0;0;0];
end

