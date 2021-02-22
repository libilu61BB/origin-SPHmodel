function [Rho_P,Rho_W] = density(person_x,person_y,wall_x,wall_y,h)
%density 计算每个粒子核近似后的密度
%   person_x 各行人粒子的x坐标
%   person_y 各行人粒子的y坐标
%   wall_x 各障碍粒子的x坐标
%   wall_y 各障碍粒子的y坐标
%   m_person 常数，各行人的质量（假设行人质量相等）
%   m_wall 常数，各障碍的质量（假设障碍质量相等）
%   h 常数，核半径
%   Rho_P 待输出的各行人粒子的密度
%   Rho_W 待输出的各障碍粒子的密度
%% 设置初始参数
n=length(person_x);
s=length(wall_x);
m_person=70;%设置行人粒子的质量为70kg
m_wall=500;%设置障碍粒子的质量为500kg
Rho_P=zeros(1,n);%初始化行人粒子的初始密度矩阵
Rho_W=zeros(1,s);%初始化障碍粒子的初始密度矩阵
%% 判断输入参数是否合法
if length(person_y)~=n
    error('行人的xy坐标数量不一致');
else
    if length(wall_y)~=s
        error('障碍的xy坐标数量不一致');
    end
end
%% 计算行人粒子的密度
%avg_Radius=mean(Radius);
%Rho_0=m_person*(4/(pi*h^8))*(h^2-4*avg_Radius^2)^3;%临界密度
for i=1:n
    for j=1:n
        %if i==j
            %continue;
        %end
        r_2=(person_x(i)-person_x(j))^2+(person_y(i)-person_y(j))^2;%计算两行人粒子之间的距离平方
        if r_2<=h^2 %只计算核半径范围内其他行人粒子对粒子i密度的影响
            Rho_P(i)=Rho_P(i)+m_person*(4/(pi*h^8))*(h^2-r_2)^3;
        end
    end
    %if Rho_P(i)<1
        %Rho_P(i)=1;
    %end
    d=zeros(1,s);
    for j=1:s
        d(j)=(person_x(i)-wall_x(j))^2+(person_y(i)-wall_y(j))^2;%计算行人与障碍粒子的距离平方
    end
    [r,u]=min(d);
    if r<=h^2 %只计算核半径范围内最近障碍粒子对粒子i密度的影响
        Rho_P(i)=Rho_P(i)+m_wall*(4/(pi*h^8))*(h^2-d(u))^3;
    end
end
%% 计算障碍粒子的密度
for i=1:s
    %for j=1:s
        %if i==j
            %continue;
        %end
        %r_2=(wall_x(i)-wall_x(j))^2+(wall_y(i)-wall_y(j))^2;%计算两障碍粒子之间的距离平方
        %if r_2<=0.01 %障碍粒子密度不受其他障碍粒子影响
    Rho_W(i)=Rho_W(i)+m_wall*(4/(pi*h^2));
        %end
    %end
    d=zeros(1,n);
    for j=1:n
        d(j)=(wall_x(i)-person_x(j))^2+(wall_y(i)-person_y(j))^2;%计算障碍粒子与行人粒子的距离平方
        if d(j)<=h^2 %只计算核半径范围内行人粒子对粒子i密度的影响
            Rho_W(i)=Rho_W(i)+m_person*(4/(pi*h^8))*(h^2-d(j))^3;
        end
    end
end
