% Author: Philipp Jost
% Purpose: U1 - drawellipse

% takes #of edges, center m, radius r, and a Matrix T to perform
% a linear Transformation of space which affects the circle.
% example
% ellipse(20,[0 0],5,[
% 1 1
% 3 1])
function ellipse(edges,m,r,T)

% use r for both m1 and m2
phi = linspace(0,2*pi,edges);
x= m(1) + r * cos(phi);
y= m(2) + r * sin(phi);

% matrix with all vectors to all points of the ellipse
A = [
    x
    y
];

TA = T * A;



plot(TA(1,:),TA(2,:));
legend(strcat('ellipse with center ',num2str(m)));
grid on
axis equal
end