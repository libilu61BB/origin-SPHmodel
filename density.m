function [Rho_P,Rho_W] = density(person_x,person_y,wall_x,wall_y,h)
%density ����ÿ�����Ӻ˽��ƺ���ܶ�
%   person_x ���������ӵ�x����
%   person_y ���������ӵ�y����
%   wall_x ���ϰ����ӵ�x����
%   wall_y ���ϰ����ӵ�y����
%   m_person �����������˵���������������������ȣ�
%   m_wall ���������ϰ��������������ϰ�������ȣ�
%   h �������˰뾶
%   Rho_P ������ĸ��������ӵ��ܶ�
%   Rho_W ������ĸ��ϰ����ӵ��ܶ�
%% ���ó�ʼ����
n=length(person_x);
s=length(wall_x);
m_person=70;%�����������ӵ�����Ϊ70kg
m_wall=500;%�����ϰ����ӵ�����Ϊ500kg
Rho_P=zeros(1,n);%��ʼ���������ӵĳ�ʼ�ܶȾ���
Rho_W=zeros(1,s);%��ʼ���ϰ����ӵĳ�ʼ�ܶȾ���
%% �ж���������Ƿ�Ϸ�
if length(person_y)~=n
    error('���˵�xy����������һ��');
else
    if length(wall_y)~=s
        error('�ϰ���xy����������һ��');
    end
end
%% �����������ӵ��ܶ�
%avg_Radius=mean(Radius);
%Rho_0=m_person*(4/(pi*h^8))*(h^2-4*avg_Radius^2)^3;%�ٽ��ܶ�
for i=1:n
    for j=1:n
        %if i==j
            %continue;
        %end
        r_2=(person_x(i)-person_x(j))^2+(person_y(i)-person_y(j))^2;%��������������֮��ľ���ƽ��
        if r_2<=h^2 %ֻ����˰뾶��Χ�������������Ӷ�����i�ܶȵ�Ӱ��
            Rho_P(i)=Rho_P(i)+m_person*(4/(pi*h^8))*(h^2-r_2)^3;
        end
    end
    %if Rho_P(i)<1
        %Rho_P(i)=1;
    %end
    d=zeros(1,s);
    for j=1:s
        d(j)=(person_x(i)-wall_x(j))^2+(person_y(i)-wall_y(j))^2;%�����������ϰ����ӵľ���ƽ��
    end
    [r,u]=min(d);
    if r<=h^2 %ֻ����˰뾶��Χ������ϰ����Ӷ�����i�ܶȵ�Ӱ��
        Rho_P(i)=Rho_P(i)+m_wall*(4/(pi*h^8))*(h^2-d(u))^3;
    end
end
%% �����ϰ����ӵ��ܶ�
for i=1:s
    %for j=1:s
        %if i==j
            %continue;
        %end
        %r_2=(wall_x(i)-wall_x(j))^2+(wall_y(i)-wall_y(j))^2;%�������ϰ�����֮��ľ���ƽ��
        %if r_2<=0.01 %�ϰ������ܶȲ��������ϰ�����Ӱ��
    Rho_W(i)=Rho_W(i)+m_wall*(4/(pi*h^2));
        %end
    %end
    d=zeros(1,n);
    for j=1:n
        d(j)=(wall_x(i)-person_x(j))^2+(wall_y(i)-person_y(j))^2;%�����ϰ��������������ӵľ���ƽ��
        if d(j)<=h^2 %ֻ����˰뾶��Χ���������Ӷ�����i�ܶȵ�Ӱ��
            Rho_W(i)=Rho_W(i)+m_person*(4/(pi*h^8))*(h^2-d(j))^3;
        end
    end
end
