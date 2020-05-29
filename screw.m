%Function to calculate screw directions for each one of shoulder, elbow and
%wrist articulations

%Input parameters:
%l = array with limb lengths: [arm, forearm and hand]
%rcg = array with proximal center of mass distances: [arm, forearm, hand]
%Alpha = array with the displaciments (rad) for each of the 7 DoF:
%   Alpha(1) = Shoulder flexion/extension
%   Alpha(2) = Shoulder abduction/adduction
%   Alpha(3) = Shoulder's internal rotation
%   Alpha(4) = Elbow's flexion/extension
%   Alpha(5) = Elbow's supination/pronation
%   Alpha(6) = Wrist's flexion/extension
%   Alpha(7) = Wrist's abduction/adduction
%r = distance of the weight's center of mass to the hand. It should be
%zero if no weight is being hold by the hand

%Output parameters
%s = screw directions [s1, s2, s3, s4, s5, s5x4, s6, s7, sx6x7]
%s0 = screw's coordinate positions
%   s0(:,1) = shoulder (s0a)
%   s0(:,2) = elbow (s0b)
%   s0(:,3) = wrist (s0c)
%   s0(:,4) = basketball's center of mass (ball joint) (s0d)
%   s0(:,5) = prismatic z virtual joint (s0e)
%   s0(:,6) = prismatic y virtual joint (s0f)
%   s0(:,7) = prismatic x virtual joint (s0g)
%   s0(:,8) = upper arm's center of mass
%   s0(:,9) = forearm's center of mass
%   s0(:,10) = wrist's center of mass


function [s, s0] = screw(l, rcg, Alpha,r)

%Basic screw directions
su = [1; 0; 0];
sv = [0; 1; 0];
sw = [0; 0; 1];

%Computes the direction vectors: s
s = zeros(3,9);
s(:,2) = rodrigues_rot(su,sw,Alpha(1));
s(:,3) = rodrigues_rot(sv,sw,Alpha(1));
s(:,3) = rodrigues_rot(s(:,3),s(:,2),Alpha(2));
s(:,1) = cross(s(:,2),s(:,3))/norm(cross(s(:,2),s(:,3)));
s(:,4) = rodrigues_rot(s(:,1),s(:,3),Alpha(3));

%s(:,1) = rodrigues_rot(s(:,1),s(:,3),Alpha(3));
%s(:,2) = rodrigues_rot(s(:,2),s(:,3),Alpha(3));

s(:,5) = rodrigues_rot(s(:,3),s(:,4),Alpha(4));
s(:,6) = cross(s(:,5),s(:,4))/norm(cross(s(:,5),s(:,4)));
s(:,7) = rodrigues_rot(s(:,4),s(:,5),Alpha(5));
s(:,9) = rodrigues_rot(s(:,5),s(:,7),Alpha(6));
s(:,8) = cross(s(:,9),s(:,7))/norm(cross(s(:,9),s(:,7)));
s(:,9) = rodrigues_rot(s(:,9),s(:,8),Alpha(7));
s(:,7) = rodrigues_rot(s(:,7),s(:,8),Alpha(7));

%Computes the position vector: s0
s0 = zeros(3,10);
s0(:,1) = [0; 0; 0];
s0(:,2) = s0(:,1) - l(1).*s(:,3);
s0(:,3) = s0(:,2) - l(2).*s(:,5);
s0(:,4) = s0(:,3) - l(3).*s(:,9) + r.*s(:,8);
s0(:,5) = [s0(1,4); s0(2,1); s0(3,4)];
s0(:,6) = [s0(1,4); s0(2,1); s0(3,1)];
s0(:,7) = [s0(1,1); s0(2,1); s0(3,1)];
s0(:,8) = s0(:,1)-(l(1)*rcg(1)).*s(:,3);
s0(:,9) = s0(:,2)-(l(2)*rcg(2)).*s(:,5);
s0(:,10) = s0(:,3)-(l(3)*rcg(3)).*s(:,9);

end

