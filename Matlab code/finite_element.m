% driver function for getting Laplacian matrix
% also mass matrix but I haven't finish that part of my code
% and eventually an approximation for u in poisson eqn

% for my matlab to work properly :/
 clc, close, clear all;

% [p,t] = icosphere(n), larger n -> more triangles on the mesh
[p,t] = icosphere(3); % code for testing FEM on a sphere (no bd)
fd=@(p) sqrt(sum(p.^2,2))-1;
% [pp,t]=distmesh2d(fd,@huniform,0.1,[-1,-1;1,1],[]); % code for testing FEM on a circle (with bd)
% p_0 = zeros(size(pp, 1), 1);
% p = [pp(:,:) p_0(:,:)];
 
% need number of vertices, meaning I want to find the largest # in t
index_max = max(t, [], 'all');
 
% construct Laplacian matrix and mass matrix
L = zeros(index_max);
M = zeros(index_max);
for v = 1:index_max
    [sumCotan, nodes_list_e] = getCotan(p, t, v);
    [area] = getArea(p, t, v);
    for n = 1:size(nodes_list_e, 1)
        L(v, nodes_list_e(n)) = -0.5*sum(sumCotan(n,:));
        M(v, nodes_list_e(n)) = area(n);
    end
end
L = L-diag(sum(L,2));
M = M+diag(sum(M,2));

% testing
theta = [];
b = [];
a = [];
a_check = [];
for i = 1:size(p, 1)
    [~, theta(i), ~] = cart2sph(p(i, 1), p(i, 2), p(i, 3));
%     b(i) = 2*p(i,1)*(p(i,3)^3)+6*p(i,1)*(p(i,2)^2)*p(i,3);
%     a_check(i) = p(i,1)*(p(i,2)^2)*(p(i,3)^3);
%     b(i) = -sqrt(pi/3)*cos(theta(i));
%     a_check(i) = 1/2*sqrt(pi/3)*cos(theta(i));
    b(i) = -1.5*sqrt(5/pi)*(3*(cos(theta(i))^2)-1);
    a_check(i) = 1/4*sqrt(5/pi)*(3*(cos(theta(i))^2)-1);
end

% err_inf = abs(L*(M*a_check')-M*b')/norm(M*b');
err = L*(M*a_check')-M*b';
err_L2 = sqrt(dot(err, M*err))

% La = Mb

% a = M*b.'\L
% diff = a-a_check

%subplot(2,1,1)
%trisurf(t,p(:,1),p(:,2),0*p(:,1),err,'edgecolor','k','facecolor','interp');
%view(2),axis([-1 1 -1 1]),axis equal, colorbar;

% subplot(2,1,2); 
% trisurf(t,p(:,1),p(:,2),0*p(:,1),a,'edgecolor','k','facecolor','interp');
% view(2),axis([-1 1 -1 1]),axis equal;