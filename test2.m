clear;
%% �����ϰ������ꡢ��������ͳ�������
%15m��15m�����οռ䣬���ڿ��3m
% wall_x1=(15:-0.1:0);wall_y1=zeros(1,length(wall_x1));
% wall_y2=(0:0.1:15);wall_x2=zeros(1,length(wall_y2));
% wall_x3=(0:0.1:15);wall_y3=15*ones(1,length(wall_x3));
% wall_y4=(15:-0.1:7.75);wall_x4=15*ones(1,length(wall_y4));
% wall_y5=(6.25:-0.1:0);wall_x5=15*ones(1,length(wall_y5));
% wall_x=[wall_x5 wall_x1 wall_x2 wall_x3 wall_x4];
% wall_y=[wall_y5 wall_y1 wall_y2 wall_y3 wall_y4];
% % �ڿռ���������ɵ�����ģ������
% % person_x=0.1+14.8*rand(1,100);
% % person_y=0.1+14.8*rand(1,100);
% load personInSquare.mat
% exit_x=16;%����x����
% end_x = 15;%�������
% exit_y=7;%����y����
%% ��ʼ��2X100mͨ��������
wall_x1 = (-100:0.1:100);
wall_y1 = zeros(1, length(wall_x1));
wall_x2 = (-100:0.1:100);
wall_y2 = 2 * ones(1, length(wall_x2));
wall_x = [wall_x1 wall_x2];
wall_y = [wall_y1 wall_y2];
load 2X100wall.mat
load personIni.mat
exit_x=101;
exit_y=1;
end_x=100;


%% �������꣬����ͼ��
n=length(person_x);
s=length(wall_x);
h1=5;%�����ܶȺ��ų���ʱʹ�õĺ˰뾶
Radius=0.3*ones(1,n);%�������˵İ뾶��Ϊ0.3m
m_person=70;%���˵�����
m_wall=500;%�ϰ��������
va=1.2; %�������˵������ٶ�ֵ
vb=1.8; %�������˵������ٶ�ֵ
u=2; %ճ�ȣ����ڼ�������֮��Ħ���������ļ��ٶ�
vx=zeros(1,n);%�����ٶ���x�����ϵķ�������ʼʱ��Ϊ0
vy=zeros(1,n);%�����ٶ���y�����ϵķ�������ʼʱ��Ϊ0
T=200; %ģ����ʱ��
sum_escape=0;%ͳ������ɢ������
P=1;%��Ϥ����·�ߵ����˱���
P_f=1;%���ڳ̶�
dt=0.03;
for t=0:dt:T
    %% �����ų����ͼ�ѹ�������ļ��ٶ�
    [Rho_person,Rho_wall]=density(person_x,person_y,wall_x,wall_y,h1);%���ú���density����tʱ�̵��ܶ�
    [ar_x,ar_y]=a_repul(person_x,person_y,wall_x,wall_y,Rho_person ,Rho_wall,Radius,h1);%���ú���a_repul�����ų��������ļ��ٶ�
    [ae_x,ae_y]=a_extru(person_x,person_y,wall_x,wall_y,Rho_person,Rho_wall,Radius,h1);%���ú���a_extru���㼷ѹ�������ļ��ٶ�
    %% ����Ħ���������ļ��ٶ�
    av_x=zeros(1,n);%��ʼ��Ħ����������x������ٶ�
    av_y=zeros(1,n);%��ʼ��Ħ����������y������ٶ�
    for i=1:n
        for j=1:n %��������������֮���Ħ���������ļ��ٶ�
            if j==i
                continue;
            end
            r=sqrt((person_x(i)-person_x(j))^2+(person_y(i)-person_y(j))^2); %����������֮��ľ���
            if r<=(Radius(i)+Radius(j))
                av_x(i)=av_x(i)+u*m_person*((vx(j)-vx(i))/(Rho_person(i)*Rho_person(j)))*(1800*(Radius(i)+Radius(j)-r))/(45*pi*(Radius(i)+Radius(j))^5);
                av_y(i)=av_y(i)+u*m_person*((vy(j)-vy(i))/(Rho_person(i)*Rho_person(j)))*(1800*(Radius(i)+Radius(j)-r))/(45*pi*(Radius(i)+Radius(j))^5);
            end
        end
        d=zeros(1,s);%��ʼ����������i����ϰ����ӵľ���
        for j=1:s
            d(j)=sqrt((person_x(i)-wall_x(j))^2+(person_y(i)-wall_y(j))^2);
        end
        [r,m]=min(d); %rΪd����Сֵ��mΪd����Сֵ������
        if r<=Radius(i) %�������ϰ�֮���Ħ���������ļ��ٶ�
            av_x(i)=av_x(i)+u*m_wall*((0-vx(i))/(Rho_person(i)*Rho_wall(m)))*(1800*(Radius(i)-r))/(45*pi*(Radius(i))^5);
            av_y(i)=av_y(i)+u*m_wall*((0-vy(i))/(Rho_person(i)*Rho_wall(m)))*(1800*(Radius(i)-r))/(45*pi*(Radius(i))^5);
        end
    end
    %% �������������������ļ��ٶ�
    e_x=zeros(1,n); %��ʼ������������x����
    e_y=zeros(1,n); %��ʼ������������y����
    am_x=zeros(1,n); %��ʼ��������������x������ٶ�
    am_y=zeros(1,n); %��ʼ��������������y������ٶ�
    for i=1:(P*n) %������Ϥ����·�ߵ����˵��˶�����
        r=sqrt((person_x(i)-exit_x)^2+(person_y(i)-exit_y)^2); %��Ϥ����·�ߵ��������������֮��ľ���
        e_x(i)=(exit_x-person_x(i))/r;
        e_y(i)=(exit_y-person_y(i))/r;
    end
    for i=(P*n+1):n %Ϊ����Ϥ����·�ߵ����˲���һ��������˶�����
        e_x(i)=-1+2*rand;
        e_y(i)=-1+2*rand;
        r=sqrt(e_x(i)^2+e_y(i)^2);%����������ģ
        e_x(i)=e_x(i)/r;%����������ת��Ϊ��λ��������
        e_y(i)=e_y(i)/r;%����������ת��Ϊ��λ��������
    end
    for i=(P*n+1):n %���㲻��Ϥ����·�ߵ������ڴ�����ΪӰ���µ��˶�����
        e_x(i)=(1-P_f)*e_x(i);
        e_y(i)=(1-P_f)*e_y(i);
        for j=1:n
            if j==i
                continue;
            end
            r_2=(person_x(i)-person_x(j))^2+(person_y(i)-person_y(j))^2; %���������Ӿ����ƽ��
            if r_2<=h1^2
                e_x(i)=e_x(i)+P_f*(m_person*e_x(j)/Rho_person(j))*(4/(pi*h1^8))*(h1^2-r_2)^3;
                e_y(i)=e_y(i)+P_f*(m_person*e_y(j)/Rho_person(j))*(4/(pi*h1^8))*(h1^2-r_2)^3;
            end
        end
        r=sqrt(e_x(i)^2+e_y(i)^2);%����������ģ
        e_x(i)=e_x(i)/r;%����������ת��Ϊ��λ��������
        e_y(i)=e_y(i)/r;%����������ת��Ϊ��λ��������
        r=sqrt((person_x(i)-exit_x)^2+(person_y(i)-exit_y)^2);
        if r<=5 %������Ϥ��ɢ·�ߵ����˽ӽ�����ʱ��������Ϊ��Ϥ��ɢ·�ߵ�����
            if person_x(i)>=exit_x
                e_x(i)=(exit_x-person_x(i))/r;
                e_y(i)=(exit_y-person_y(i))/r;
            end
        end
    end
    for i=1:n
%         am_x(i)=(v0*e_x(i)-vx(i))/dt;
%         am_y(i)=(v0*e_y(i)-vy(i))/dt;        
        am_x(i)=min((va*e_x(i)-vx(i))/dt,20);
        am_y(i)=min((va*e_y(i)-vy(i))/dt,20);
    end
    %% �������˵�λ��
    ax = am_x+ar_x+ae_x+av_x;%1��n�У�tʱ�̸���������x����ĺϼ��ٶ�
    ay = am_y+ar_y+ae_y+av_y;%1��n�У�tʱ�̸���������y����ĺϼ��ٶ�
    vx = vx+ax*dt; %������һʱ�̵�x�����ٶ�
    vy = vy+ay*dt; %������һʱ�̵�y�����ٶ�
    V = sqrt(vx.^2+vy.^2);
    V1 = V(1:floor(n/2));%ǰһ������˵ĺ��ٶȣ�floorΪ����ȡ��
    V2 = V(floor(n/2)+1:n);%��һ�����˵ĺ��ٶ�
    index_V1 = find(V1>va);%�ҳ����ٵĵ������ˣ�������������
    index_V2 = find(V2>vb);%�ҳ����ٵĸ������ˣ�������������
    
    vx(index_V1) = vx(index_V1).*va./V1(index_V1);%����ǰһ������Ϊ��������
    vy(index_V1) = vy(index_V1).*va./V1(index_V1);
    vx(index_V2+floor(n/2)) = vx(index_V2+floor(n/2)).*vb./V2(index_V2);%�����һ������Ϊ��������
    vy(index_V2+floor(n/2)) = vy(index_V2+floor(n/2)).*vb./V2(index_V2); 
%     for i=1:n %����һʱ�̵��ٶȴ���v0��������С��v0
%         r = sqrt(vx(i)^2+vy(i)^2);
%         if r>=va
%             vx(i) = vx(i)*va/r;
%             vy(i) = vy(i)*va/r;
%         end
%     end
%     for i=1:n
%         if person_x(i)>15
%             vx(i) = va;
%             vy(i) = 0;
%         end
%     end
    person_x = person_x+vx*dt; %����x�����λ��
    person_y = person_y+vy*dt; %����y�����λ��
    
    %%
    [minP2P, minP2W] = minDistance(person_x,person_y,wall_x,wall_y);
%     if minP2P < 0.5
%         fprintf("person min < 0.5");
%     end
%     if minP2W < 0.3
%         fprintf("wall min < 0.3");
%     end
    
    %%
    %for i=1:n
        %if person_x(i)>exit_x %����ͨ�����ں���x�����뿪
            %vx(i)=v0;
            %vy(i)=0;
            %continue;
        %end       
    %end 
    %% ����ͼ��
    for i=1:n %ͳ����������
        if person_x(i)>(end_x)
            sum_escape=sum_escape+1;
            person_x(i)=nan;
            person_y(i)=nan;
        end
    end
    
%     % 15X15��ͼ
%     plot(person_x,person_y,'.',wall_x,wall_y,'LineWidth',2);
%     hold on;
%     plot(person_x,person_y,'.', 'MarkerSize', 10)
%     axis([-1 18 -1 18]);%������ʾ��Χ
%     set(gcf,'position',[0,0,1000,1000]);
    
    % 2X100��ͼ
    plot(wall_x1,wall_y1,'LineWidth',1,'Color','k');
    hold on;
    plot(wall_x2,wall_y2,'LineWidth',1,'Color','k');
    hold on;
    plot(person_x,person_y,'.', 'MarkerSize', 10)
    axis([-1 101 -2 4]);%������ʾ��Χ
%     set(gcf,'position',[0,0,1000,1000]);

    str_time=sprintf('��ɢʱ�䣺%.2f',t);
    str_escape=sprintf('����������%.0f',sum_escape);
    text(5.5,-0.5,str_time);
    text(5.5,-1.3,str_escape);
    axis on;
    hold off;
    pause(0.001);
        if sum_escape>=n
            break;
        end
end
        
        
        
    
                
                