function [ae_x,ae_y] = a_extru(person_x,person_y,wall_x,wall_y,Rho_person,Rho_wall,Radius,disP2P,disP2W,m_person,m_wall)
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
n=length(person_x);
ae_x=zeros(1,n);
ae_y=zeros(1,n);
ae_p2pMax = 80; %行人之间挤压力产生的加速度的最大值
h2 = 2*Radius;
h3 = Radius;

%% 计算各粒子的压强
Rho_p2p=2;
Rho_p2w=2;
Pe_person=P_extru(Rho_person,Rho_p2p); %调用函数P_extru计算行人粒子的压强
Pe_wall=P_extru(Rho_wall,Rho_p2w); %调用函数P_extru计算障碍粒子的压强

%% 计算行人粒子的加速度分量
for i=1:n
    % ↓↓计算行人粒子之间的挤压加速度↓↓
    disp2p = disP2P(i,1:n); %粒子i与其他粒子的距离
    disp2p(i) = nan; %不考虑粒子与自己的距离
    index = find(disp2p<h2); %找到发生挤压的粒子的索引
    m = length(index);
    Pe_person_i = Pe_person(i)*ones(1,m); %将粒子i的压强转化为与index维度一致的矩阵
    Rho_person_i = Rho_person(i)*ones(1,m); %将粒子i的密度转化为与index维度一致的矩阵
    abs_aeP2P = m_person*(Pe_person_i./Rho_person_i.^2+Pe_person(index)./Rho_person(index).^2).*...
        (3*(10*(h2-disp2p(index)).^2)/(pi*h2^5));
    person_xi = person_x(i)*ones(1,m);
    person_yi = person_y(i)*ones(1,m);
    ae_x(i) = sum(abs_aeP2P.*(person_xi-person_x(index))./disp2p(index));
    ae_y(i) = sum(abs_aeP2P.*(person_yi-person_y(index))./disp2p(index));
    % ↓↓若行人粒子之间的挤压力加速度过大，则将其缩小↓↓
    ae = sqrt(ae_x(i)^2+ae_y(i)^2);
    if ae>ae_p2pMax
        ae_x(i) = ae_x(i)*ae_p2pMax/ae;
        ae_y(i) = ae_y(i)*ae_p2pMax/ae;
    end   
end
% ↓↓计算行人与障碍之间的挤压加速度↓↓
[disp2w,ind_w] = min(disP2W(:,1:n)); %disP2W按列计算最小值，并返回每列最小值的行号，相当于障碍物粒子的索引
ind_p = find(disp2w<h3); %计算与障碍物发生碰撞的行人粒子的索引
if ~isempty(ind_p) %如果有粒子与障碍物碰撞，则计算人与障碍物之间的挤压加速度
    abs_aeP2W = m_wall*(Pe_person(ind_p)./Rho_person(ind_p).^2+Pe_wall(ind_w(ind_p))./Rho_wall(ind_w(ind_p)).^2).*...
        (3*(10*(h3-disp2w(ind_p)).^2)/(pi*h3^5));
    ae_x(ind_p)=ae_x(ind_p)+abs_aeP2W.*(person_x(ind_p)-wall_x(ind_w(ind_p)))./disp2w(ind_p); %将加速度分解到x轴
    ae_y(ind_p)=ae_y(ind_p)+abs_aeP2W.*(person_y(ind_p)-wall_y(ind_w(ind_p)))./disp2w(ind_p); %将加速度分解到y轴
end