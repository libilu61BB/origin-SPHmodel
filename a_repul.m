function [ar_x,ar_y] = a_repul(person_x,person_y,wall_x,wall_y,Rho_person,Rho_wall,Radius,h1)
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
m_person=70; %���˵�����
m_wall=500; %�ϰ�������
n=length(person_x);
s=length(wall_x);
ar_x=zeros(1,n);
ar_y=zeros(1,n);
%% �ж���������Ƿ�Ϸ�
if length(person_y)~=n
    error('���˵�xy����������һ��');
else
    if length(wall_y)~=s
        error('�ϰ���xy����������һ��');
    end
end
%% ��������Ӻ�ѹǿ
avg_Radius=mean(Radius);
h2 = 2*avg_Radius;
% Rho_p2p=m_person*(4/(pi*h1^8))*(h1^2-4*avg_Radius^2)^3+m_person*(4/(pi*h1^2)); %������֮����ٽ��ܶ�
% % Rho_p2p=m_person*(4/(pi*h1^8))*(h1^2-4*avg_Radius^2)^3; %������֮����ٽ��ܶ�
% Rho_p2w=m_person*(4/(pi*h1^8))*(h1^2-avg_Radius^2)^3+m_wall*(4/(pi*h1^2)); %�����ϰ�֮����ٽ��ܶ�
% % Rho_p2w=m_wall*(4/(pi*h1^8))*(h1^2-avg_Radius^2)^3; %�����ϰ�֮����ٽ��ܶ�

% test
Rho_p2p=2;
Rho_p2w=2;
Pr_person=P_repul(Rho_person,Rho_p2p); %���ú���P_repul�����������ӵ�ѹǿ
Pr_wall=P_repul(Rho_wall,Rho_p2w); %���ú���P_repul�����ϰ����ӵ�ѹǿ
%% �����������ӵļ��ٶȷ���   
for i=1:n
    %-----------------��������������֮����ų��������ļ��ٶ�-----------------
    for j=1:n
        if j==i
            continue;
        end
        r=sqrt((person_x(i)-person_x(j))^2+(person_y(i)-person_y(j))^2); %����������֮��ľ���
        if r<=h1
            abs_ar=m_person*(Pr_person(i)/Rho_person(i)^2+Pr_person(j)/Rho_person(j)^2)*(3*(10*(h1-r)^2)/(pi*h1^5));
            %���ٶȵ�ģ��
            ar_x(i)=ar_x(i)+abs_ar*(person_x(i)-person_x(j))/r; %�����ٶȷֽ⵽x��
            ar_y(i)=ar_y(i)+abs_ar*(person_y(i)-person_y(j))/r; %�����ٶȷֽ⵽y��
        end
    end
    %-----------------�����������ϰ�֮����ų��������ļ��ٶ�-----------------
    d=zeros(1,s);
    for j=1:s
        d(j)=sqrt((person_x(i)-wall_x(j))^2+(person_y(i)-wall_y(j))^2);
    end
    [r,u]=min(d); %rΪd����Сֵ��uΪd����Сֵ������
    if r<=h2
        abs_ar=m_wall*(Pr_person(i)/Rho_person(i)^2+Pr_wall(u)/Rho_wall(u)^2)*(3*(10*(h2-r)^2)/(pi*h2^5));
        %���ٶȵ�ģ��
        ar_x(i)=ar_x(i)+abs_ar*(person_x(i)-wall_x(u))/r; %�����ٶȷֽ⵽x��
        ar_y(i)=ar_y(i)+abs_ar*(person_y(i)-wall_y(u))/r; %�����ٶȷֽ⵽y��
    end
end
end
