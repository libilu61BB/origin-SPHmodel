function [ae_x,ae_y] = a_extru(person_x,person_y,wall_x,wall_y,Rho_person,Rho_wall,Radius,disP2P,disP2W,m_person,m_wall)
%a_extru �������Ӽ��໥�Ӵ�ʱ��ѹ�������ļ��ٶȵ�xy����
%   person_x �������ӵ�x����
%   person_y �������ӵ�y����
%   wall_x �ϰ����ӵ�x����
%   wall_y �ϰ����ӵ�y����
%   Radius ���������ӵİ뾶
%   h1 ���������������ܶȺ��ų���ʱʹ�õĺ˰뾶
%   h2 ���������㼷ѹ�������ļ��ٶ�ʱ���õĺ˰뾶
%   ae_x ���ٶ���x�����ϵķ���
%   ae_y ���ٶ���y�����ϵķ���
%   Pe_person �Ӵ��������ļ�ѹ���Ը��������Ӳ�����ѹǿ
%   Pe_wall �Ӵ��������ļ�ѹ���Ը��ϰ����Ӳ�����ѹǿ
%   Rho_person ���������ӵ��ܶ�
%   Rho_wall ���ϰ����ӵ��ܶ�
%% ���ó�ʼ����
n=length(person_x);
ae_x=zeros(1,n);
ae_y=zeros(1,n);
ae_p2pMax = 80; %����֮�伷ѹ�������ļ��ٶȵ����ֵ
h2 = 2*Radius;
h3 = Radius;

%% ��������ӵ�ѹǿ
Rho_p2p=2;
Rho_p2w=2;
Pe_person=P_extru(Rho_person,Rho_p2p); %���ú���P_extru�����������ӵ�ѹǿ
Pe_wall=P_extru(Rho_wall,Rho_p2w); %���ú���P_extru�����ϰ����ӵ�ѹǿ

%% �����������ӵļ��ٶȷ���
for i=1:n
    % ����������������֮��ļ�ѹ���ٶȡ���
    disp2p = disP2P(i,1:n); %����i���������ӵľ���
    disp2p(i) = nan; %�������������Լ��ľ���
    index = find(disp2p<h2); %�ҵ�������ѹ�����ӵ�����
    m = length(index);
    Pe_person_i = Pe_person(i)*ones(1,m); %������i��ѹǿת��Ϊ��indexά��һ�µľ���
    Rho_person_i = Rho_person(i)*ones(1,m); %������i���ܶ�ת��Ϊ��indexά��һ�µľ���
    abs_aeP2P = m_person*(Pe_person_i./Rho_person_i.^2+Pe_person(index)./Rho_person(index).^2).*...
        (3*(10*(h2-disp2p(index)).^2)/(pi*h2^5));
    person_xi = person_x(i)*ones(1,m);
    person_yi = person_y(i)*ones(1,m);
    ae_x(i) = sum(abs_aeP2P.*(person_xi-person_x(index))./disp2p(index));
    ae_y(i) = sum(abs_aeP2P.*(person_yi-person_y(index))./disp2p(index));
    % ��������������֮��ļ�ѹ�����ٶȹ���������С����
    ae = sqrt(ae_x(i)^2+ae_y(i)^2);
    if ae>ae_p2pMax
        ae_x(i) = ae_x(i)*ae_p2pMax/ae;
        ae_y(i) = ae_y(i)*ae_p2pMax/ae;
    end   
end
% ���������������ϰ�֮��ļ�ѹ���ٶȡ���
[disp2w,ind_w] = min(disP2W(:,1:n)); %disP2W���м�����Сֵ��������ÿ����Сֵ���кţ��൱���ϰ������ӵ�����
ind_p = find(disp2w<h3); %�������ϰ��﷢����ײ���������ӵ�����
if ~isempty(ind_p) %������������ϰ�����ײ������������ϰ���֮��ļ�ѹ���ٶ�
    abs_aeP2W = m_wall*(Pe_person(ind_p)./Rho_person(ind_p).^2+Pe_wall(ind_w(ind_p))./Rho_wall(ind_w(ind_p)).^2).*...
        (3*(10*(h3-disp2w(ind_p)).^2)/(pi*h3^5));
    ae_x(ind_p)=ae_x(ind_p)+abs_aeP2W.*(person_x(ind_p)-wall_x(ind_w(ind_p)))./disp2w(ind_p); %�����ٶȷֽ⵽x��
    ae_y(ind_p)=ae_y(ind_p)+abs_aeP2W.*(person_y(ind_p)-wall_y(ind_w(ind_p)))./disp2w(ind_p); %�����ٶȷֽ⵽y��
end