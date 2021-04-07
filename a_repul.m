function [ar_x,ar_y] = a_repul(person_x,person_y,wall_x,wall_y,Rho_person,Rho_wall,Radius,disP2P,disP2W,m_person,m_wall,h1)
%a_repul �����ų����ڸ����������ϲ����ļ��ٶ�
%   person_x ���������ӵ�x����
%   person_y ���������ӵ�y����
%   wall_x ���ϰ����ӵ�x����
%   wall_y ���ϰ����ӵ�y����
%   h1 ���������������ܶȺ��ų���ʱʹ�õĺ˰뾶
%   ar_x ���ٶ���x�����ϵķ���
%   ar_y ���ٶ���y�����ϵķ���
%   Pr_person �ų����Ը��������Ӳ�����ѹǿ
%   Pr_wall �ų����Ը��ϰ����Ӳ�����ѹǿ
%   Rho_person ���������ӵ��ܶ�
%   Rho_wall ���ϰ����ӵ��ܶ�

%% ���ó�ʼ����
n=length(person_x);
s=length(wall_x);
ar_x=zeros(1,n);
ar_y=zeros(1,n);
h2 = 2*Radius;

%% ��������ӵ�ѹǿ
Rho_p2p = 2;
Rho_p2w = 2;
Pr_person = P_repul(Rho_person,Rho_p2p); %���ú���P_repul�����������ӵ�ѹǿ
Pr_wall = P_repul(Rho_wall,Rho_p2w); %���ú���P_repul�����ϰ����ӵ�ѹǿ

%% �����������ӵļ��ٶȷ���   
for i=1:n
    % ����������������֮����ų���ٶȡ���
    disp2p = disP2P(i,1:n); %����i���������ӵľ���
    disp2p(i) = nan; %�������������Լ��ľ���
    index = find(disp2p<h1); %�ҵ�������ѹ�����ӵ�����
    m = length(index);
    Pr_person_i = Pr_person(i)*ones(1,m); %������i��ѹǿת��Ϊ��indexά��һ�µľ���
    Rho_person_i = Rho_person(i)*ones(1,m); %������i���ܶ�ת��Ϊ��indexά��һ�µľ���
    abs_arP2P = m_person*(Pr_person_i./Rho_person_i.^2+Pr_person(index)./Rho_person(index).^2).*...
        (3*(10*(h1-disp2p(index)).^2)/(pi*h1^5));
    person_xi = person_x(i)*ones(1,m);
    person_yi = person_y(i)*ones(1,m);
    ar_x(i) = sum(abs_arP2P.*(person_xi-person_x(index))./disp2p(index));
    ar_y(i) = sum(abs_arP2P.*(person_yi-person_y(index))./disp2p(index));
end
% ���������������ϰ�֮����ų���ٶȡ���
[disp2w,ind_w] = min(disP2W(:,1:n)); %disP2W���м�����Сֵ��������ÿ����Сֵ���кţ��൱���ϰ������ӵ�����
ind_p = find(disp2w<h2); %�������ϰ��﷢���ų���������ӵ�����
if ~isempty(ind_p) %������������ϰ�������ų⣬����������ϰ���֮����ų���ٶ�
    abs_arP2W = m_wall*(Pr_person(ind_p)./Rho_person(ind_p).^2+Pr_wall(ind_w(ind_p))./Rho_wall(ind_w(ind_p)).^2).*...
        (3*(10*(h2-disp2w(ind_p)).^2)/(pi*h2^5));
    ar_x(ind_p) = ar_x(ind_p)+abs_arP2W.*(person_x(ind_p)-wall_x(ind_w(ind_p)))./disp2w(ind_p); %�����ٶȷֽ⵽x��
    ar_y(ind_p) = ar_y(ind_p)+abs_arP2W.*(person_y(ind_p)-wall_y(ind_w(ind_p)))./disp2w(ind_p); %�����ٶȷֽ⵽y��
end
