function [ar_x,ar_y] = a_repul(person_x,person_y,wall_x,wall_y,Rho_person,Rho_wall,Radius,disP2P,disP2W,m_person,m_wall,h1)
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
n=length(person_x);
s=length(wall_x);
ar_x=zeros(1,n);
ar_y=zeros(1,n);
h2 = 2*Radius;

%% 计算各粒子的压强
Rho_p2p = 2;
Rho_p2w = 2;
Pr_person = P_repul(Rho_person,Rho_p2p); %调用函数P_repul计算行人粒子的压强
Pr_wall = P_repul(Rho_wall,Rho_p2w); %调用函数P_repul计算障碍粒子的压强

%% 计算行人粒子的加速度分量   
for i=1:n
    % ↓↓计算行人粒子之间的排斥加速度↓↓
    disp2p = disP2P(i,1:n); %粒子i与其他粒子的距离
    disp2p(i) = nan; %不考虑粒子与自己的距离
    index = find(disp2p<h1); %找到发生挤压的粒子的索引
    m = length(index);
    Pr_person_i = Pr_person(i)*ones(1,m); %将粒子i的压强转化为与index维度一致的矩阵
    Rho_person_i = Rho_person(i)*ones(1,m); %将粒子i的密度转化为与index维度一致的矩阵
    abs_arP2P = m_person*(Pr_person_i./Rho_person_i.^2+Pr_person(index)./Rho_person(index).^2).*...
        (3*(10*(h1-disp2p(index)).^2)/(pi*h1^5));
    person_xi = person_x(i)*ones(1,m);
    person_yi = person_y(i)*ones(1,m);
    ar_x(i) = sum(abs_arP2P.*(person_xi-person_x(index))./disp2p(index));
    ar_y(i) = sum(abs_arP2P.*(person_yi-person_y(index))./disp2p(index));
end
% ↓↓计算行人与障碍之间的排斥加速度↓↓
[disp2w,ind_w] = min(disP2W(:,1:n)); %disP2W按列计算最小值，并返回每列最小值的行号，相当于障碍物粒子的索引
ind_p = find(disp2w<h2); %计算与障碍物发生排斥的行人粒子的索引
if ~isempty(ind_p) %如果有粒子与障碍物产生排斥，则计算人与障碍物之间的排斥加速度
    abs_arP2W = m_wall*(Pr_person(ind_p)./Rho_person(ind_p).^2+Pr_wall(ind_w(ind_p))./Rho_wall(ind_w(ind_p)).^2).*...
        (3*(10*(h2-disp2w(ind_p)).^2)/(pi*h2^5));
    ar_x(ind_p) = ar_x(ind_p)+abs_arP2W.*(person_x(ind_p)-wall_x(ind_w(ind_p)))./disp2w(ind_p); %将加速度分解到x轴
    ar_y(ind_p) = ar_y(ind_p)+abs_arP2W.*(person_y(ind_p)-wall_y(ind_w(ind_p)))./disp2w(ind_p); %将加速度分解到y轴
end
