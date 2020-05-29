%Function to calculate the Bezier curve, given the points, the derivatives
%and the number of discretization points for each section of the curve

%Input parameters:
%r = array will all points in which the Bezier curve shall pass
%rl0 = derivative at the initial point
%rln = derivative at the final point
%N = number of points in which each section of the curve is discretized

%Output parameters:
%y = points of the Bezier curve
%yl = first derivative of the points of the Bezier curve
%yll = second derivative of the points of the Bezier curve

function [y, yl, yll] = bezier(r, rl0, rln, N)

%Verifies the number of points of the Bezier
n = length(r);

%Creates the matrix where the Bezier curve will be assembled: M*V(i) = C(i)
M = zeros(2*(n-1),2*(n-1));

%Insets two additional lines to the system of equations regarding first 
%and final derivatives
M(2*(n-1)-1,1) = 1;
M(2*(n-1),2*(n-1)) = 1;
C(2*(n-1)-1,1) = (3*r(1)+rl0)/3;
C(2*(n-1),1) = (3*r(n)-rln)/3;

%Computes all control points for the Bezier
V = M\C;

%Divides each section of the Bezier into N points
u = [0:(1/N):1];

%Determines the points of the Bezier 
%Each column of y, yl, yll corresponds to a section of the curve
for i=1:n-1
    
    %Computes the matrix H used for calculation of Bezier curves
    H = [1 0 0 0;-3 3 0 0;3 -6 3 0;-1 3 -3 1];
    
    for j = 1:(length(u)-1)
        
        %Calculates the Bezier points y
        U = [1, u(j), u(j)^2, u(j)^3];
        Y = [r(i); V(2*i-1); V(2*i); r(i+1)];
        y(j,i) = U*H*Y;
        
        %Calculates Bezier first derivative yl
        UL = 3.*[(1-u(j))^2, 2*u(j)*(1-u(j)), u(j)^2];
        YL = [V(2*i-1)-r(i); V(2*i)-V(2*i-1); r(i+1)-V(2*i)];
        yl(j,i) = UL*YL;
        
        %Calculates Bezier second derivative yll
        ULL = 6.*[1-u(j), u(j)];
        YLL = [Y(3) - 2*Y(2) + Y(1); Y(2) - 2*Y(3) + Y(4)];
        yll(j,i) = ULL*YLL;
        
    end
    
end

%Converts the matrices y, yll and yll in only one array
y = y(:);
yl = yl(:);
yll = yll(:);

%Adds one last element of the curve that was not calculated through the for
%loop above
j = length(u);
U = [1 u(j) u(j)^2 u(j)^3];
Y = [r(i) V(2*i-1) V(2*i) r(i+1)]';
UL = 3.*[(1-u(j))^2, 2*u(j)*(1-u(j)), u(j)^2];
YL = [V(2*i-1)-r(i); V(2*i)-V(2*i-1); r(i+1)-V(2*i)];
ULL = 6.*[1-u(j), u(j)];
YLL = [Y(3) - 2*Y(2) + Y(1); Y(2) - 2*Y(3) + Y(4)];
y(length(y)+1) = U*H*Y;
yl(length(yl)+1) = UL*YL;
yll(length(yll)+1) = ULL*YLL;


end

