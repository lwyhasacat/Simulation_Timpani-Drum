function [pts,tris]=cylinder(cnt)

%cnt = 100;
uvmesh = 1/cnt;
topmesh = uvmesh*2;
u = uvmesh:uvmesh:1;
v = 0:uvmesh:1;
v2 = topmesh:topmesh:1-topmesh;
theta = 2*pi * u;
x1 = cos(theta);
%x2 = v2'*cos(theta);
x2 = reshape((v2'*cos(theta)).', 1, []);
y1 = sin(theta);
y2 = reshape((v2'*sin(theta)).', 1, []);
z = -2 * v + 1;
psizeside = cnt*(cnt+1); % cnt^2 + cnt points on the side
psizetop = cnt * (cnt/2-1); % 0.5*cnt^2-cnt points on one top
psize = psizeside + 2 * psizetop + 2;  % 2cnt^2 - cnt + 2
pts = zeros(psize, 3);
for i = 1:psizeside
    pts(i, 1) = x1(rem(i-1, cnt)+1);
    pts(i, 2) = y1(rem(i-1, cnt)+1);
    pts(i, 3) = z(fix((i-1)/cnt)+1); % division with remainder
end
pts(psizeside+1, :) = [0,0,1]; % top center point
for i = psizeside+2: psizeside+psizetop+1
    pts(i, 1) = x2(i-psizeside-1);
    pts(i, 2) = y2(i-psizeside-1);
    pts(i, 3) = 1;
end
pts(psizeside+psizetop+2, :) = [0, 0, -1]; % bottom center point
for i = psizeside+psizetop+3: psizeside+2*psizetop+2
    pts(i, 1) = x2(i-psizeside-psizetop-2);
    pts(i, 2) = y2(i-psizeside-psizetop-2);
    pts(i, 3) = -1;
end
%disp(pts);
tsizeside = 2*cnt*cnt;
tsizetop = cnt*cnt-cnt;%0.5*cnt*cnt;
tsize = tsizeside + 2 * tsizetop;
tris = zeros(tsize, 3);
for i = 1:tsizeside
    if rem(i, 2*cnt) ~= 0 && rem(i, 2*cnt) ~= 2*cnt-1
        fd = cnt*fix(i/2/cnt);
        tris(i, 1) = fd + floor(rem(i, 2*cnt)*0.5)+1;
        tris(i, 2) = tris(i, 1) + cnt; % vertical link
        if rem(i, 2) == 1 % odd
            tris(i, 3) = tris(i,1)+1;
        else % even
            tris(i, 3) = tris(i-1,2);
        end
    elseif rem(i, 2*cnt) == 2*cnt-1
        tris(i, 1) = tris(i-1, 1);
        tris(i, 2) = tris(i, 1) + cnt; % vertical link
        tris(i, 3) = tris(i-2*cnt+2, 1);
    elseif rem(i, 2*cnt) == 0
        tris(i, 1) = tris(i-2*cnt+1, 1);
        tris(i, 2) = tris(i, 1) + cnt; % vertical link
        tris(i, 3) = tris(i-1, 2);
    end
end
for i = tsizeside+1 : tsizeside+cnt % top inner loop
    tris(i, 1) = psizeside+1;
    tris(i, 2) = psizeside+2 + i - tsizeside - 1;
    if i ~= tsizeside+cnt
        tris(i, 3) = tris(i, 2) + 1;
    else
        tris(i, 3) = psizeside+2;
    end
end
for i = tsizeside + cnt + 1: tsizeside + tsizetop-2*cnt % top middle loops
    j = i - tsizeside - cnt;
    if rem(j, 2 * cnt) ~= 0 && rem(j, 2 * cnt) ~= 2 * cnt - 1
        fd = cnt * fix(j / 2 / cnt);
        tris(i, 1) = psizeside+1 +fd + floor(rem(j, 2 * cnt) * 0.5) + 1;
        tris(i, 2) = tris(i, 1) + cnt; % vertical link
        if rem(j, 2) == 1 % odd
            tris(i, 3) = tris(i, 1)+1;
        else % even
            tris(i, 3) = tris(i-1, 2);
        end
    elseif rem(j, 2*cnt) == 2*cnt-1
        tris(i, 1) = tris(i-1, 1);
        tris(i, 2) = tris(i, 1) + cnt; % vertical link
        tris(i, 3) = tris(i-2*cnt+2, 1);
    elseif rem(j, 2*cnt) == 0
        tris(i, 1) = tris(i-2*cnt+1, 1);
        tris(i, 2) = tris(i, 1) + cnt; % vertical link
        tris(i, 3) = tris(i-1, 2);
    end
end

for i = tsizeside + tsizetop - 2*cnt + 1: tsizeside + tsizetop % top outer loop
    lpcnt = 2*cnt;
    j = i - tsizeside - tsizetop + lpcnt;
    if rem(j, lpcnt) ~= 0 && rem(j, lpcnt) ~= 2 * cnt - 1
        tris(i, 1) = psizeside+1 + (cnt/2-2)* cnt + floor(rem(j, lpcnt) * 0.5) + 1;
        tris(i, 2) = tris(i, 1) - tris(tsizeside+tsizetop-lpcnt+1, 1)+1; % vertical link
        if rem(j, 2) == 1 % odd
            tris(i, 3) = tris(i, 1)+1;
        else % even
            tris(i, 3) = tris(i-1, 2);
        end
    elseif rem(j, lpcnt) == lpcnt-1
        tris(i, 1) = tris(i-1, 1);
        tris(i, 2) = tris(i, 1) - tris(tsizeside+tsizetop-lpcnt+1, 1)+1; % vertical link
        tris(i, 3) = tris(i-lpcnt+2, 1);
    elseif rem(j, lpcnt) == 0
        tris(i, 1) = tris(i-lpcnt+1, 1);
        tris(i, 2) = tris(i, 1) - tris(tsizeside+tsizetop-lpcnt+1, 1)+1; % vertical link
        tris(i, 3) = tris(i-1, 2);
    end
end


for i = tsizeside+tsizetop+1 : tsizeside+tsizetop+cnt % bottom inner loop
    tris(i, 1) = psizeside+psizetop+1+1;
    tris(i, 2) = psizeside+psizetop+1+1 + i - tsizeside - tsizetop;
    if i ~= tsizeside+tsizetop+cnt
        tris(i, 3) = tris(i, 2) + 1;
    else
        tris(i, 3) = psizeside+psizetop+3;
    end
end
for i = tsizeside +tsizetop+ cnt + 1: tsizeside + 2*tsizetop-2*cnt % top middle loops
    j = i - tsizeside -tsizetop- cnt;
    if rem(j, 2 * cnt) ~= 0 && rem(j, 2 * cnt) ~= 2 * cnt - 1
        fd = cnt * fix(j / 2 / cnt);
        tris(i, 1) = psizeside+psizetop+1+1 +fd + floor(rem(j, 2 * cnt) * 0.5) + 1;
        tris(i, 2) = tris(i, 1) + cnt; % vertical link
        if rem(j, 2) == 1 % odd
            tris(i, 3) = tris(i, 1)+1;
        else % even
            tris(i, 3) = tris(i-1, 2);
        end
    elseif rem(j, 2*cnt) == 2*cnt-1
        tris(i, 1) = tris(i-1, 1);
        tris(i, 2) = tris(i, 1) + cnt; % vertical link
        tris(i, 3) = tris(i-2*cnt+2, 1);
    elseif rem(j, 2*cnt) == 0
        tris(i, 1) = tris(i-2*cnt+1, 1);
        tris(i, 2) = tris(i, 1) + cnt; % vertical link
        tris(i, 3) = tris(i-1, 2);
    end
end

for i = tsizeside + 2*tsizetop - 2*cnt + 1: tsizeside + 2*tsizetop % top outer loop
    lpcnt = 2*cnt;
    j = i - tsizeside - 2*tsizetop + lpcnt;
    if rem(j, lpcnt) ~= 0 && rem(j, lpcnt) ~= 2 * cnt - 1
        tris(i, 1) = psizeside+1 + psizetop+1+(cnt/2-2)* cnt + floor(rem(j, lpcnt) * 0.5) + 1;
        tris(i, 2) = tris(i, 1) - tris(tsizeside+2*tsizetop-lpcnt+1, 1)+cnt*cnt+1; % vertical link
        if rem(j, 2) == 1 % odd
            tris(i, 3) = tris(i, 1)+1;
        else % even
            tris(i, 3) = tris(i-1, 2);
        end
    elseif rem(j, lpcnt) == lpcnt-1
        tris(i, 1) = tris(i-1, 1);
        tris(i, 2) = tris(i, 1) - tris(tsizeside+2*tsizetop-lpcnt+1, 1)+cnt*cnt+1; % vertical link
        tris(i, 3) = tris(i-lpcnt+2, 1);
    elseif rem(j, lpcnt) == 0
        tris(i, 1) = tris(i-lpcnt+1, 1);
        tris(i, 2) = tris(i, 1) - tris(tsizeside+2*tsizetop-lpcnt+1, 1)+cnt*cnt+1; % vertical link
        tris(i, 3) = tris(i-1, 2);
    end
end