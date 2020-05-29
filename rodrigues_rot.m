% Rodrigues' rotation formula
% Direction is determined by the right-hand (screw) rule.
%
% Inputs:
%    v     - Array of three dimensional vectors to rotate. Array can be 
%            composed of N rows of 3D row vectors or N columns of 3D
%            column vectors. If v is 3x3 array, it is assumed that it is 3
%            rows of 3 3D row vectors
%    k     - Rotation axis (does not need to be unit vector)
%    theta - Rotation angle in radians; positive according to right-hand
%           (screw) rule
%
%   Note: k and individual 3D vectors in v array must be same orientation.
%           
% Outputs:
%    v_rot - Array of rotated vectors.


function v_rot = rodrigues_rot(v,k,theta)
    
    %Checks if there is any error in the input parameters
    [m,n] = size(v);
    if (m ~= 3 && n ~= 3)
        error('input vector is/are not three dimensional'), end
    if (size(v) ~= size(k)) 
        error('rotation vector v and axis k have different dimensions'),end
    
    %Normalizes the rotation axis
    k = k/sqrt(k(1)^2 + k(2)^2 + k(3)^2);
    
    %Verifies the number vectors in the array
    No = numel(v)/3;
    
    %Initializes the rotated vector array
    v_rot = v;
    
    %If the input v represents row vectors (m ~= 3 && n == 3)
    if ( n == 3 )
        %Initializes the cross product k and v with right dim
        crosskv = v(1,:); 
        for i = 1:No
            crosskv(1) = k(2)*v(i,3) - k(3)*v(i,2);
            crosskv(2) = k(3)*v(i,1) - k(1)*v(i,3); 
            crosskv(3) = k(1)*v(i,2) - k(2)*v(i,1);
            v_rot(i,:) = cos(theta)*v(i,:) + (crosskv)*sin(theta)...
                            + k*(dot(k,v(i,:)))*(1 - cos(theta));
        end
    %If the input v represents column vectors (n ~= 3 && m == 3)    
    else
        %Initializes the cross product k and v with right dim
        crosskv = v(:,1);
        for i = 1:No
            crosskv(1) = k(2)*v(3,i) - k(3)*v(2,i);
            crosskv(2) = k(3)*v(1,i) - k(1)*v(3,i); 
            crosskv(3) = k(1)*v(2,i) - k(2)*v(1,i);
            v_rot(:,i) = cos(theta)*v(:,i) + (crosskv)*sin(theta)...
                            + k*(dot(k,v(:,i)))*(1 - cos(theta));
        end
    end
end

