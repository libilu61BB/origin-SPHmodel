function [Rho_P,Rho_W] = density(n, m_person, m_wall, h, disP2P, disP2W)
%density 计算每个粒子核近似后的密度
%   n 行人粒子的数量
%   s 障碍粒子的数量
%   m_person 常数，各行人的质量（假设行人质量相等）
%   m_wall 常数，各障碍的质量（假设障碍质量相等）
%   h 常数，核半径
%   disP2P 行人之间的距离
%   disP2W 行人与障碍之间的距离
%   Rho_P 待输出的各行人粒子的密度
%   Rho_W 待输出的各障碍粒子的密度

% % 设置初始参数
% Rho_P=zeros(1,n);%初始化行人粒子的初始密度矩阵
% Rho_W=zeros(1,s);%初始化障碍粒子的初始密度矩阵

%% 计算行人粒子的密度
% ↓↓行人与行人↓↓
RhoP2P_temp = m_person*(4/(pi*h^8))*(h^2-disP2P(1:n,1:n).^2).^3; %计算行人ij之间的密度贡献
RhoP2P_temp(RhoP2P_temp<0) = 0; %把负密度（说明不在核半径内）重置为0
RhoP2P = sum(RhoP2P_temp); %按列求和，得到1×n的密度矩阵，每一列代表对应粒子的P2P密度
% ↓↓行人与障碍物↓↓
disP2W_min = min(disP2W(:,1:n));
RhoP2W_temp = m_wall*(4/(pi*h^8))*(h^2-disP2W_min(1).^2).^3; %按列求disP2W的最小值，得到每个行人与障碍物的最小距离，计算距离最近的障碍物对行人粒子的密度贡献
RhoP2W_temp(RhoP2W_temp<0) = 0; %把负密度（说明不在核半径内）重置为0
Rho_P2W = RhoP2W_temp;
% ↓↓密度求和↓↓
Rho_P = RhoP2P+Rho_P2W;

%% 计算障碍粒子的密度
RhoW2P_temp = m_person*(4/(pi*h^8))*(h^2-disP2W(:,1:n).^2).^3; %计算行人粒子对障碍物密度的贡献
RhoW2P_temp(RhoW2P_temp<0) = 0;
Rho_W = m_wall*(4/(pi*h^2))+sum(RhoW2P_temp,2)';