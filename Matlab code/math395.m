% driver function for solving Laplacian problems on a bunny
clc, close, clear all;

% save video
video = true; % if want to save a movie, then let video = true
if(video)
    writerObj = VideoWriter('drum_sim_damping.mp4','MPEG-4');
    writerObj.FrameRate = 20;
    open(writerObj);
end

% get p, t matrices
% [p,t]=cylinder(20);
[p,t]=cylinder(20);

p = p./max(p(:));
dh = norm(p(t(1,1),:)-p(t(1,2),:));

% wave speed
lamb = 3.5676;
freq = 130.81;
c = freq*lamb;

dt = 2*10^(-5);

% get Laplacian and Mass matrices
[L,M] = getLM(t,p);

h = trisurf(t, p(:,1), p(:,2), p(:,3));
daspect([1 1 1])
axis([-2 2 -2 2 -2 2])
drawnow

v = 0*p;
theta = atan2(sqrt(p(:,1).^2 + p(:,2).^2),p(:,3));
pt = 0.5*repmat(exp(-10*(theta).^2),1,3);
Tf = 1.0;
NtimeSteps = 50000;

% plot every 'nplot' time steps
nplot = 100;
S=zeros(1, NtimeSteps);
for i = 1:NtimeSteps
    k0 = 0.000002;
    for k = 1:3
        v(:,k) = v(:,k) - dt*c^2*(M\(L*pt(:,k) + k0*v(:, k))); % damping
        % v(:,k) = v(:,k) - dt*c^2*(M\(L*pt(:,k)));
        pt(:,k) = pt(:,k) + dt*v(:,k);
    end
    S(1, i) = norm(pt(1,:));
    if mod(i,nplot) == 0
        pl = p+pt;
        h = trisurf(t, pl(:,1), pl(:,2), pl(:,3));
        daspect([1 1 1])
        axis([-2 2 -2 2 -2 2])
        title(['time = ' num2str((i*dt))])
        drawnow
        if video
            frame = getframe(gcf);
            writeVideo(writerObj,frame);
        end
    end
end
figure(2);
plot(1:NtimeSteps, S);
if video
    close(writerObj);
end
soundsc(S,NtimeSteps)

filename = 'sim_C3_betterMesh.wav';
audiowrite(filename,S,NtimeSteps);


function [L,M] = getLM(t,p)
index_max = max(t, [], 'all');

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
end