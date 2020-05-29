%Function to plot arm posture and the screws of shoulder, elbow and wrist
%
%Input parameters:
%l = array with upper arm, arm and hand lengths
%rcg = center of mass positions for the upper arm, arm and hand
%Alpha = rotation for each one of the shoulder, elbow ans wrist DoF's
%   Alpha(1) = Shoulder flexion/extension
%   Alpha(2) = Shoulder abduction/adduction
%   Alpha(3) = Shoulder's internal rotation
%   Alpha(4) = Elbow's flexion/extension
%   Alpha(5) = Elbow's supination/pronation
%   Alpha(6) = Wrist's flexion/extension
%   Alpha(7) = Wrist's abduction/adduction
%r = basketball's radius

function [] = plotposture(l,rcg, Alpha,r)

%Determines screw directions and articulation positions
[s,s0] = screw(l,rcg, Alpha,r);

figure
hold on;
%Plots articulation coordinates
plot3([s0(1,1) s0(1,2)],[s0(2,1) s0(2,2)],[s0(3,1) s0(3,2)],'k');
plot3([s0(1,2) s0(1,3)],[s0(2,2) s0(2,3)],[s0(3,2) s0(3,3)],'r');
plot3([s0(1,3) s0(1,4)],[s0(2,3) s0(2,4)],[s0(3,3) s0(3,4)],'b');
%Plots center of masses positions
plot3(s0(1,8),s0(2,8),s0(3,8),'*k')
plot3(s0(1,9),s0(2,9),s0(3,9),'*r')
plot3(s0(1,10),s0(2,10),s0(3,10),'*b')
%Configures the graph
xlabel('Eixo X');
ylabel('Eixo Y');
zlabel('Eixo Z');
grid on;
%axis([-1 1 -1 1 -1 1])
%axis([0.12 0.2 -0.02 0.06 0 0.2])
axis([-0.3 0.8 -0.4 0.7 -0.4 0.7])

end

