function [Rho_P,Rho_W] = density(n, m_person, m_wall, h, disP2P, disP2W)
%density ����ÿ�����Ӻ˽��ƺ���ܶ�
%   n �������ӵ�����
%   s �ϰ����ӵ�����
%   m_person �����������˵���������������������ȣ�
%   m_wall ���������ϰ��������������ϰ�������ȣ�
%   h �������˰뾶
%   disP2P ����֮��ľ���
%   disP2W �������ϰ�֮��ľ���
%   Rho_P ������ĸ��������ӵ��ܶ�
%   Rho_W ������ĸ��ϰ����ӵ��ܶ�

% % ���ó�ʼ����
% Rho_P=zeros(1,n);%��ʼ���������ӵĳ�ʼ�ܶȾ���
% Rho_W=zeros(1,s);%��ʼ���ϰ����ӵĳ�ʼ�ܶȾ���

%% �����������ӵ��ܶ�
% �������������ˡ���
RhoP2P_temp = m_person*(4/(pi*h^8))*(h^2-disP2P(1:n,1:n).^2).^3; %��������ij֮����ܶȹ���
RhoP2P_temp(RhoP2P_temp<0) = 0; %�Ѹ��ܶȣ�˵�����ں˰뾶�ڣ�����Ϊ0
RhoP2P = sum(RhoP2P_temp); %������ͣ��õ�1��n���ܶȾ���ÿһ�д����Ӧ���ӵ�P2P�ܶ�
% �����������ϰ������
disP2W_min = min(disP2W(:,1:n));
RhoP2W_temp = m_wall*(4/(pi*h^8))*(h^2-disP2W_min(1).^2).^3; %������disP2W����Сֵ���õ�ÿ���������ϰ������С���룬�������������ϰ�����������ӵ��ܶȹ���
RhoP2W_temp(RhoP2W_temp<0) = 0; %�Ѹ��ܶȣ�˵�����ں˰뾶�ڣ�����Ϊ0
Rho_P2W = RhoP2W_temp;
% �����ܶ���͡���
Rho_P = RhoP2P+Rho_P2W;

%% �����ϰ����ӵ��ܶ�
RhoW2P_temp = m_person*(4/(pi*h^8))*(h^2-disP2W(:,1:n).^2).^3; %�����������Ӷ��ϰ����ܶȵĹ���
RhoW2P_temp(RhoW2P_temp<0) = 0;
Rho_W = m_wall*(4/(pi*h^2))+sum(RhoW2P_temp,2)';