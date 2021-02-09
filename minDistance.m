function [minP2P, minP2W] = minDistance(person_x,person_y,wall_x,wall_y)
%minDistance 计算各个粒子之间的最小距离
%   person_x 各行人粒子的x坐标
%   person_y 各行人粒子的y坐标
%   wall_x 各障碍粒子的x坐标
%   wall_y 各障碍粒子的y坐标
%   minP2P 行人粒子之间的最小距离
%   minP2W 行人与障碍粒子之间的最小距离
n = length(person_x);
s = length(wall_x);
minP2P = 10000;
minP2W = 10000;
for i=1:n
    for j =1:n
        if j==i
            continue;
        end
        r=sqrt((person_x(i)-person_x(j))^2+(person_y(i)-person_y(j))^2); %两行人粒子之间的距离
        if r < minP2P
            minP2P = r;
        end
    end
    for j = 1:s
        d = sqrt((person_x(i)-wall_x(j))^2+(person_y(i)-wall_y(j))^2);
        if d < minP2W
            minP2W = d;
        end
    end
end