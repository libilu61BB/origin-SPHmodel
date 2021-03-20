function [person_l_num, person_l_x, person_l_y, person_r_num, person_r_x, person_r_y] = personInitialization(density_l,density_r,l_width,l_height,r_width,r_height)
% person initialization
% 根据不同的密度完成行人初始化
% density_l 通道左侧需要达到的行人密度
% density_r 通道右侧需要达到的行人密度
% l_width 左侧初始化区域宽度
% l_height 左侧初始化区域高度
% r_width 右侧初始化区域宽度
% r_height 右侧初始化区域高度

person_l_x = [];
person_l_y = [];
person_r_x = [];
person_r_y = [];
%% 我要开始初始化左边惹！
person_l_num = ceil(density_l * l_width * l_height);
while length(person_l_x) < person_l_num
    flag=1;
    temp_x = 0.25 + (l_width-0.5) * rand(1,1);
    temp_y = 0.25 + (l_height-0.5) * rand(1,1);
    for i=1:length(person_l_x)
        r = sqrt((temp_x-person_l_x(i))^2 + (temp_y-person_l_y(i))^2);
        if r < 0.44
            flag = 0;
            break;
        end
    end
    if flag == 1
        person_l_x = [person_l_x temp_x];
        person_l_y = [person_l_y temp_y];
    end
end
% % 啊画一下左边初始化的行人粒子示意图康康~
% plot(person_l_x, person_l_y, '.');
% axis([0 l_width 0 l_height]);%设置显示范围
% set(gcf,'position',[200,200,100*l_width,100*l_height]);

%% 我要开始初始化右边惹！
person_r_num = ceil(density_r * r_width * r_height);
while length(person_r_x) < person_r_num
    flag=1;
    temp_x = 0.25 + (r_width-0.5) * rand(1,1);
    temp_y = 0.25 + (r_height-0.5) * rand(1,1);
    for i=1:length(person_r_x)
        r = sqrt((temp_x-person_r_x(i))^2 + (temp_y-person_r_y(i))^2);
        if r < 0.4
            flag = 0;
            break;
        end
    end
    if flag == 1
        person_r_x = [person_r_x temp_x];
        person_r_y = [person_r_y temp_y];
    end
end
% % 啊画一下右边初始化的行人粒子示意图康康~
% plot(person_r_x, person_r_y, '.r');
% axis([0 r_width 0 r_height]);%设置显示范围
% set(gcf,'position',[200,200,100*l_width,100*l_height]);
