function [ar_x,ar_y] = a_repul(person_x,person_y,wall_x,wall_y,Rho_person,Rho_wall,Radius,h1)
%a_repul 计算排斥力在各行人粒子上产生的加速度
%   person_x 各行人粒子的x坐标
%   person_y 各行人粒子的y坐标
%   wall_x 各障碍粒子的x坐标
%   wall_y 各障碍粒子的y坐标
%   h1 常数，计算粒子密度和排斥力时使用的核半径
%   ar_x 加速度在x方向上的分量
%   ar_y 加速度在y方向上的分量
%   Pr_person 排斥力对各行人粒子产生的压强
%   Pr_wall 排斥力对各障碍粒子产生的压强
%   Rho_person 各行人粒子的密度
%   Rho_wall 各障碍粒子的密度
%% 设置初始参数
m_person=70; %行人的质量
m_wall=500; %障碍的质量
n=length(person_x);
s=length(wall_x);
ar_x=zeros(1,n);
ar_y=zeros(1,n);
%% 判断输入参数是否合法
if length(person_y)~=n
    error('行人的xy坐标数量不一致');
else
    if length(wall_y)~=s
        error('障碍的xy坐标数量不一致');
    end
end
%% 计算各粒子和压强
avg_Radius=mean(Radius);
h2 = 2*avg_Radius;
% Rho_p2p=m_person*(4/(pi*h1^8))*(h1^2-4*avg_Radius^2)^3+m_person*(4/(pi*h1^2)); %人与人之间的临界密度
% % Rho_p2p=m_person*(4/(pi*h1^8))*(h1^2-4*avg_Radius^2)^3; %人与人之间的临界密度
% Rho_p2w=m_person*(4/(pi*h1^8))*(h1^2-avg_Radius^2)^3+m_wall*(4/(pi*h1^2)); %人与障碍之间的临界密度
% % Rho_p2w=m_wall*(4/(pi*h1^8))*(h1^2-avg_Radius^2)^3; %人与障碍之间的临界密度

% test
Rho_p2p=2;
Rho_p2w=2;
Pr_person=P_repul(Rho_person,Rho_p2p); %调用函数P_repul计算行人粒子的压强
Pr_wall=P_repul(Rho_wall,Rho_p2w); %调用函数P_repul计算障碍粒子的压强
%% 计算行人粒子的加速度分量   
for i=1:n
    %-----------------计算行人与行人之间的排斥力产生的加速度-----------------
    for j=1:n
        if j==i
            continue;
        end
        r=sqrt((person_x(i)-person_x(j))^2+(person_y(i)-person_y(j))^2); %两行人粒子之间的距离
        if r<=h1
            abs_ar=m_person*(Pr_person(i)/Rho_person(i)^2+Pr_person(j)/Rho_person(j)^2)*(3*(10*(h1-r)^2)/(pi*h1^5));
            %加速度的模长
            ar_x(i)=ar_x(i)+abs_ar*(person_x(i)-person_x(j))/r; %将加速度分解到x轴
            ar_y(i)=ar_y(i)+abs_ar*(person_y(i)-person_y(j))/r; %将加速度分解到y轴
        end
    end
    %-----------------计算行人与障碍之间的排斥力产生的加速度-----------------
    d=zeros(1,s);
    for j=1:s
        d(j)=sqrt((person_x(i)-wall_x(j))^2+(person_y(i)-wall_y(j))^2);
    end
    [r,u]=min(d); %r为d中最小值，u为d中最小值的索引
    if r<=h2
        abs_ar=m_wall*(Pr_person(i)/Rho_person(i)^2+Pr_wall(u)/Rho_wall(u)^2)*(3*(10*(h2-r)^2)/(pi*h2^5));
        %加速度的模长
        ar_x(i)=ar_x(i)+abs_ar*(person_x(i)-wall_x(u))/r; %将加速度分解到x轴
        ar_y(i)=ar_y(i)+abs_ar*(person_y(i)-wall_y(u))/r; %将加速度分解到y轴
    end
end
end
