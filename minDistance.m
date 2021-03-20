function [minP2P, minP2W] = minDistance(person_x,person_y,wall_x,wall_y)
%minDistance �����������֮�����С����
%   person_x ���������ӵ�x����
%   person_y ���������ӵ�y����
%   wall_x ���ϰ����ӵ�x����
%   wall_y ���ϰ����ӵ�y����
%   minP2P ��������֮�����С����
%   minP2W �������ϰ�����֮�����С����
n = length(person_x);
s = length(wall_x);
minP2P = 10000;
minP2W = 10000;
for i=1:n
    for j =1:n
        if j==i
            continue;
        end
        r=sqrt((person_x(i)-person_x(j))^2+(person_y(i)-person_y(j))^2); %����������֮��ľ���
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