function [ae_x,ae_y] = a_extru(person_x,person_y,wall_x,wall_y,Rho_person,Rho_wall,Radius,h1)
%a_extru 计算粒子间相互接触时挤压力产生的加速度的xy分量
%   person_x 行人粒子的x坐标
%   person_y 行人粒子的y坐标
%   wall_x 障碍粒子的x坐标
%   wall_y 障碍粒子的y坐标
%   Radius 各行人粒子的半径
%   h1 常数，计算粒子密度和排斥力时使用的核半径
%   h2 常数，计算挤压力产生的加速度时所用的核半径
%   ae_x 加速度在x方向上的分量
%   ae_y 加速度在y方向上的分量
%   Pe_person 接触所产生的挤压力对各行人粒子产生的压强
%   Pe_wall 接触所产生的挤压力对各障碍粒子产生的压强
%   Rho_person 各行人粒子的密度
%   Rho_wall 各障碍粒子的密度
%% 设置初始参数
m_person=70; %行人的质量
m_wall=500; %障碍的质量
n=length(person_x);
s=length(wall_x);
ae_x=zeros(1,n);
ae_y=zeros(1,n);
ae_p2pMax = 80; %行人之间挤压力产生的加速度的最大值
%% 判断输入参数是否合法
if length(person_y)~=n
    error('行人的xy坐标数量不一致');
else
    if length(wall_y)~=s
        error('障碍的xy坐标数量不一致');
    end
end
%% 计算各粒子和压强
% avg_Radius=mean(Radius);
% Rho_p2p=m_person*(4/(pi*h1^8))*(h1^2-4*avg_Radius^2)^3+m_person*(4/(pi*h1^2)); %人与人之间的临界密度
% Rho_p2w=m_person*(4/(pi*h1^8))*(h1^2-avg_Radius^2)^3+m_wall*(4/(pi*h1^2)); %人与障碍之间的临界密度
% Rho_p2p=m_person*(4/(pi*h1^8))*(h1^2-4*avg_Radius^2)^3; %人与人之间的临界密度
% Rho_p2w=m_person*(4/(pi*h1^8))*(h1^2-avg_Radius^2)^3; %人与障碍之间的临界密度

% test
Rho_p2p=2;
Rho_p2w=2;
Pe_person=P_extru(Rho_person,Rho_p2p); %调用函数P_extru计算行人粒子的压强
Pe_wall=P_extru(Rho_wall,Rho_p2w); %调用函数P_extru计算障碍粒子的压强
%% 计算行人粒子的加速度分量
for i=1:n
    %---------------行人与行人之间相互接触时挤压力产生的加速度---------------
    for j=1:n
        if j==i
            continue;
        end
        r=sqrt((person_x(i)-person_x(j))^2+(person_y(i)-person_y(j))^2); %两行人粒子之间的距离
        if r<=(Radius(i)+Radius(j)) %当两行人粒子的距离小于等于其半径之和时视为接触
            h2=Radius(i)+Radius(j);
            abs_ae=m_person*(Pe_person(i)/Rho_person(i)^2+Pe_person(j)/Rho_person(j)^2)*(3*(10*(h2-r)^2)/(pi*h2^5));
            %加速度的模长
            ae_x(i)=ae_x(i)+abs_ae*(person_x(i)-person_x(j))/r; %将加速度分解到x轴
            ae_y(i)=ae_y(i)+abs_ae*(person_y(i)-person_y(j))/r; %将加速度分解到y轴
        end
    end
    ae = sqrt(ae_x(i)^2+ae_y(i)^2);
    if ae>ae_p2pMax %若挤压力加速度过大，则将其缩小
        ae_x(i) = ae_x(i)*ae_p2pMax/ae;
        ae_y(i) = ae_y(i)*ae_p2pMax/ae;
    end
    
    %---------------行人与障碍之间相互接触时挤压力产生的加速度---------------
    d=zeros(1,s);
    for j=1:s
        d(j)=sqrt((person_x(i)-wall_x(j))^2+(person_y(i)-wall_y(j))^2);
    end
    [r,u]=min(d); %r为d中最小值，u为d中最小值的索引
    if r<=Radius(i) %当两粒子的距离小于等于其半径之和时视为接触，其中障碍粒子的半径为0
        h2=Radius(i);
        abs_ae=m_wall*(Pe_person(i)/Rho_person(i)^2+Pe_wall(u)/Rho_wall(u)^2)*(3*(10*(h2-r)^2)/(pi*h2^5));
        ae_x(i)=ae_x(i)+abs_ae*(person_x(i)-wall_x(u))/r; %将加速度分解到x轴
        ae_y(i)=ae_y(i)+abs_ae*(person_y(i)-wall_y(u))/r; %将加速度分解到y轴
    end
end
end

