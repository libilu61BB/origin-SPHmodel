function [ae_x,ae_y] = a_extru(person_x,person_y,wall_x,wall_y,Rho_person,Rho_wall,Radius,h1)
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
m_person=70; %���˵�����
m_wall=500; %�ϰ�������
n=length(person_x);
s=length(wall_x);
ae_x=zeros(1,n);
ae_y=zeros(1,n);
ae_p2pMax = 80; %����֮�伷ѹ�������ļ��ٶȵ����ֵ
%% �ж���������Ƿ�Ϸ�
if length(person_y)~=n
    error('���˵�xy����������һ��');
else
    if length(wall_y)~=s
        error('�ϰ���xy����������һ��');
    end
end
%% ��������Ӻ�ѹǿ
% avg_Radius=mean(Radius);
% Rho_p2p=m_person*(4/(pi*h1^8))*(h1^2-4*avg_Radius^2)^3+m_person*(4/(pi*h1^2)); %������֮����ٽ��ܶ�
% Rho_p2w=m_person*(4/(pi*h1^8))*(h1^2-avg_Radius^2)^3+m_wall*(4/(pi*h1^2)); %�����ϰ�֮����ٽ��ܶ�
% Rho_p2p=m_person*(4/(pi*h1^8))*(h1^2-4*avg_Radius^2)^3; %������֮����ٽ��ܶ�
% Rho_p2w=m_person*(4/(pi*h1^8))*(h1^2-avg_Radius^2)^3; %�����ϰ�֮����ٽ��ܶ�

% test
Rho_p2p=2;
Rho_p2w=2;
Pe_person=P_extru(Rho_person,Rho_p2p); %���ú���P_extru�����������ӵ�ѹǿ
Pe_wall=P_extru(Rho_wall,Rho_p2w); %���ú���P_extru�����ϰ����ӵ�ѹǿ
%% �����������ӵļ��ٶȷ���
for i=1:n
    %---------------����������֮���໥�Ӵ�ʱ��ѹ�������ļ��ٶ�---------------
    for j=1:n
        if j==i
            continue;
        end
        r=sqrt((person_x(i)-person_x(j))^2+(person_y(i)-person_y(j))^2); %����������֮��ľ���
        if r<=(Radius(i)+Radius(j)) %�����������ӵľ���С�ڵ�����뾶֮��ʱ��Ϊ�Ӵ�
            h2=Radius(i)+Radius(j);
            abs_ae=m_person*(Pe_person(i)/Rho_person(i)^2+Pe_person(j)/Rho_person(j)^2)*(3*(10*(h2-r)^2)/(pi*h2^5));
            %���ٶȵ�ģ��
            ae_x(i)=ae_x(i)+abs_ae*(person_x(i)-person_x(j))/r; %�����ٶȷֽ⵽x��
            ae_y(i)=ae_y(i)+abs_ae*(person_y(i)-person_y(j))/r; %�����ٶȷֽ⵽y��
        end
    end
    ae = sqrt(ae_x(i)^2+ae_y(i)^2);
    if ae>ae_p2pMax %����ѹ�����ٶȹ���������С
        ae_x(i) = ae_x(i)*ae_p2pMax/ae;
        ae_y(i) = ae_y(i)*ae_p2pMax/ae;
    end
    
    %---------------�������ϰ�֮���໥�Ӵ�ʱ��ѹ�������ļ��ٶ�---------------
    d=zeros(1,s);
    for j=1:s
        d(j)=sqrt((person_x(i)-wall_x(j))^2+(person_y(i)-wall_y(j))^2);
    end
    [r,u]=min(d); %rΪd����Сֵ��uΪd����Сֵ������
    if r<=Radius(i) %�������ӵľ���С�ڵ�����뾶֮��ʱ��Ϊ�Ӵ��������ϰ����ӵİ뾶Ϊ0
        h2=Radius(i);
        abs_ae=m_wall*(Pe_person(i)/Rho_person(i)^2+Pe_wall(u)/Rho_wall(u)^2)*(3*(10*(h2-r)^2)/(pi*h2^5));
        ae_x(i)=ae_x(i)+abs_ae*(person_x(i)-wall_x(u))/r; %�����ٶȷֽ⵽x��
        ae_y(i)=ae_y(i)+abs_ae*(person_y(i)-wall_y(u))/r; %�����ٶȷֽ⵽y��
    end
end
end

