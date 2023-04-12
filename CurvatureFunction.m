function [k1, k2, Gauss, Mean, dir1, dir2] = CurvatureFunction(data, query, n) %local cloud xyz data, query point xyz data, normal to the query point
coeffs = 3; %number of coefficients - form f(x,y) = 1/2(ax^2+2bxy+cy^2) as suggested by Curvature in Triangle Meshes Chapter 8
random = [1*n(1)+1, -1*n(2)-1, 4*n(1)+4]; %this transformation creates a non parallel vector to n, even where n is the zero vector
intermediate = random - n*(dot(n,random)); 
a = intermediate/sqrt(intermediate(1)^2+intermediate(2)^2+intermediate(3)^2); %find first vector defining tangent plane
b = cross(n,a); %find second vector defining the tangent plane
data = data-query; %center data at (0, 0)
T = [a', b']; %create transformation T
newxy = (T'*(data)')'; %apply transformation to centered data
u = newxy(:,1); 
v = newxy(:,2);
F = zeros(size(data(:,1)));
for subpt = 1:height(F)
    F(subpt) = dot(n,data(subpt,:)); %define new height values as n dot data
end
U = zeros(height(F), coeffs); %form coefficient matrix
U(:,1) = (u.^2)/2;
U(:,2) = u.*v;
U(:,3) = (v.^2)/2;
D = (inv(U'*U))*U'; %create least squares coefficient
X = D*F; %find coefficient matrix
a = X(1); b = X(2); c = X(3);
syms f x y
f(x,y) = (1/2)*(a*x^2 + 2*b*x*y + c*y^2); 
Shape = double(-hessian(f)); %Find the shape operator
[dir, K] = eig(Shape); %take the principal curvatures
k1 = K(1,1); dir1 = dir(:, 1); %take curvatures and directions (directions are in local CS)
k2 = K(2,2); dir2 = dir(:, 2);
K = [k1,k2]; %this section sorts principal curvatures and directions into min and max
DIR = [dir1,dir2];
[~, idx1] = max(abs(K));
[~, idx2] = min(abs(K));
k1 = K(idx1);
k2 = K(idx2);
dir1 = DIR(:,idx1);
dir2 = DIR(:,idx2);
Mean = (k1+k2)/2; %take mean curvature
Gauss = k1*k2; %take Gaussian curvature
end