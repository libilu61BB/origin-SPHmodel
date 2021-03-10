% ˫��������������Ϊģ�⣬4*50mͨ��������˫������
clear;
condition = 1; %ѡ��ģ�ⳡ��
%% ��ʼ������
wall_x1 = (-20:0.1:70);
wall_y1 = zeros(1, length(wall_x1));
wall_x2 = (-20:0.1:70);
wall_y2 = 4 * ones(1, length(wall_x2));
wall_x = [wall_x1, wall_x2];
wall_y = [wall_y1, wall_y2];
s=length(wall_x);

h1=5; %�����ܶȺ��ų���ʱʹ�õĺ˰뾶
m_person=70; %���˵�����
m_wall=500; %�ϰ��������
u=2; %ճ�ȣ����ڼ�������֮��Ħ���������ļ��ٶ�
T=200; %ģ����ʱ��
sum_escape=0; %ͳ������ɢ������
P=1; %��Ϥ����·�ߵ����˱���
P_f=1; %���ڳ̶�
dt=0.02; %ʱ�䲽��
t_person = 0.5; %�������ӵ�ʱ����

person_x = []; %���˵�x����
person_y = []; %���˵�y����
exit_x = []; %���ڵ�x����
exit_y = []; %���ڵ�y����
end_x = []; %����������
vx = []; %�����ٶ���x�����ϵķ���
vy = []; %�����ٶ���y�����ϵķ���
v0 = []; %���˵������ٶ�
Radius = []; %���˵İ뾶

%% ��֮��������������
P_r=0.5^dt;%���ٶ���֮�������ʱ��Ȩ��
A=10;%���ٶ���֮�����������

%% �ٶ�-�ܶ�ͳ������
count = 0; %�ٶ��ܶ�ͳ�Ƽ���
v_sum = 0; %ƽ���ٶȺ�
density_sum = 0; %ƽ���ܶȺ�
v_blue = [];
density_area = [];

%% ģ��ѭ��
for t=0:dt:T
    %% ���������������
    if mod(t,t_person)==0 %ÿ��0.5s�������һ��
        person_x_temp = [-0.5-5*rand 50+5*rand]; %���������������һ������
        person_y_temp = [0.3+3.4*rand 0.3+3.4*rand];
        person_x = [person_x person_x_temp];
        person_y = [person_y person_y_temp];
        exit_x = [exit_x 100 -50];
        exit_y = [exit_y 2 2];
        end_x = [end_x 50 0];
        vx = [vx 0 0];
        vy = [vy 0 0];
        temp1 = rand;
        temp2 = rand;
        v0 = [v0 1.2*(temp1>0.5)+1.8*(temp1<0.5) 1.2*(temp2>0.5)+1.8*(temp2<0.5)]; %�������˺͵��������������
        Radius = [Radius 0.3 0.3];
    end
    n=length(person_x);
    %% ��֮���������ؼ���
    al=A*rand(1,n);%���˼��ٶȵ���֮�����������Ϊ���Ӹ�˹�ֲ��������
    al_theta=2*pi*rand(1,n);%���˼��ٶ���֮����������ķ���Ϊ[0,2*pi]�ڵ������
    al_x=al.*cos(al_theta);%���˼��ٶ���x�����ϵ���֮�����������
    al_y=al.*sin(al_theta);%���˼��ٶ���y�����ϵ���֮�����������
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
            r=sqrt((person_x(i) - person_x(j))^2 + (person_y(i) - person_y(j))^2); %����������֮��ľ���
            if r<=(Radius(i)+Radius(j))
                av_x(i) = av_x(i)+u*m_person*((vx(j)-vx(i))/(Rho_person(i)*Rho_person(j)))*(1800*(Radius(i)+Radius(j)-r))/(45*pi*(Radius(i)+Radius(j))^5);
                av_y(i) = av_y(i)+u*m_person*((vy(j)-vy(i))/(Rho_person(i)*Rho_person(j)))*(1800*(Radius(i)+Radius(j)-r))/(45*pi*(Radius(i)+Radius(j))^5);
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
        r=sqrt((person_x(i)-exit_x(i))^2+(person_y(i)-exit_y(i))^2); %��Ϥ����·�ߵ��������������֮��ľ���
        e_x(i)=(exit_x(i)-person_x(i))/r;
        e_y(i)=(exit_y(i)-person_y(i))/r;
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
        r=sqrt((person_x(i)-exit_x(i))^2+(person_y(i)-exit_y(i))^2);
        if r<=5 %������Ϥ��ɢ·�ߵ����˽ӽ�����ʱ��������Ϊ��Ϥ��ɢ·�ߵ�����
            if person_x(i)>=exit_x(i)
                e_x(i)=(exit_x(i)-person_x(i))/r;
                e_y(i)=(exit_y(i)-person_y(i))/r;
            end
        end
    end
    for i=1:n
        %         am_x(i)=(v0*e_x(i)-vx(i))/dt;
        %         am_y(i)=(v0*e_y(i)-vy(i))/dt;
        temp_am_x=(v0(i)*e_x(i)-vx(i))/dt;
        if temp_am_x>=0
            am_x(i)=min((v0(i)*e_x(i)-vx(i))/dt,20);
        else
            am_x(i)=max((v0(i)*e_x(i)-vx(i))/dt,-20);
        end
        temp_am_y=(v0(i)*e_y(i)-vy(i))/dt;
        if temp_am_y>=0
            am_y(i)=min((v0(i)*e_y(i)-vy(i))/dt,20);
        else
            am_y(i)=max((v0(i)*e_y(i)-vy(i))/dt,-20);
        end
    end
    %% ������֮��������ٶ�
    yita=A*rand(1,n);%���˼��ٶȵ���֮�����������Ϊ���Ӹ�˹�ֲ��������
    theta=2*pi*rand(1,n);%���˼��ٶ���֮����������ķ���Ϊ[0,2*pi]�ڵ������
    al_x=P_r*al_x + (1-P_r)*yita.*cos(theta);%���˼��ٶ���x�����ϵ���֮�����������
    al_y=P_r*al_y + (1-P_r)*yita.*sin(theta);%���˼��ٶ���y�����ϵ���֮�����������
    %% �������˸�����Ϊ�����ļ��ٶ�
    a_graX = zeros(1,n); %��ʼ��X������ٶ�
    a_graY = zeros(1,n); %��ʼ��Y������ٶ�
    for i=1:n
        for j=1:n
            if i==j
                continue;
            end
            Dij = [person_x(j)-person_x(i),person_y(j)-person_y(i)]; %��iָ��j��λ������Dij
            Vi = [vx(i),vy(i)]; %����i���ٶ�����Vi
            Dij_abs = sqrt(sum(Dij.^2)); %����ij��ģ���൱�������ӵľ���
            if Dij_abs<=5 && sum(Dij.*Vi)>0 %������С�ڵ���2��jλ��i��ǰ��ʱ�Ž��к�������
                ei = [exit_x(i)-person_x(i),exit_y(i)-person_y(i)]/sqrt(sum([exit_x(i)-person_x(i),exit_y(i)-person_y(i)].^2)); %������iָ����ڵĵ�λ����
                Vj = [vx(j),vy(j)]; %����j���ٶ�����Vj
                Bij_3 = max(0,sum(ei.*Vj)/sqrt(sum(Vj.^2)));
                if Bij_3==0
                    continue;
                else
                    Bij_4 = min(sqrt(sum(Vj.^2))/v0(i),1);
                    Bij_5 = min(exp(1)^((Radius(i)+Radius(j))-Dij_abs),1);
                    a_graX(i) = a_graX(i)+0.4*v0(i)*Bij_3*Bij_4*Bij_5*Dij(1)/Dij_abs; %����X������ٶ�
                    a_graY(i) = a_graY(i)+0.4*v0(i)*Bij_3*Bij_4*Bij_5*Dij(2)/Dij_abs; %����Y������ٶ�
                end
            end
        end
    end
    
    %% ���㳬����Ϊ�����ļ��ٶ�
    Pl = 0; %�������ĵ÷�
    Pm = 0; %�м�����ĵ÷�
    Pr = 0; %�Ҳ�����ĵ÷�
    w_p = 0.3; %����÷�Ȩ��
    w_sa = 0.4; %ֱ��ǰ����Ȩ��
    w_rl = 0.3; %���ҳ�������õ�Ȩ��
    search_R = 5; %��������÷�ʱ�������뾶
    a_pass_abs = 20; %������Ϊ�����ļ��ٶȵĴ�С
    a_pass_x = zeros(1,n);
    a_pass_y = zeros(1,n);
    d_sa = 3;
    C_rl = 1.5;
    C_ot = 1.25;
    Kin = -1;
    h = 4;
    for i=1:n
        Vi = [vx(i),vy(i)]; %����i���ٶ�����
        Vi_abs = sqrt(sum(Vi.^2)); %����i���ٶȴ�С
        if Vi_abs==0
            continue
        end
        for j=1:n
            if i==j
                continue
            end
            Dij = [person_x(j)-person_x(i), person_y(j)-person_y(i)]; %��iָ��j��λ������Dij
            Dij_abs = sqrt(sum(Dij.^2)); %ij֮��ľ���
            if Dij_abs<=search_R
                VxD = Vi(1)*Dij(2)-Vi(2)*Dij(1); %Vi��Dij�Ĳ��
                switch (VxD>=0)
                    case 1 %VxD>=0˵��Dij��Vi����ʱ�뷽�򣬼�j��i�����
                        cos_VD = sum(Dij.*Vi)/(Dij_abs*Vi_abs);
                        if cos_VD<0.5 %�нǴ���60�㲻����
                            continue
                        else
                            if cos_VD<(sqrt(3)/2) %�н���30~60֮�䣬����Pl����
                                Vj = [vx(j),vy(j)]; %����j���ٶ�����Vj
                                Ui = [exit_x(i)-person_x(i),exit_y(i)-person_y(i)]/sqrt(sum([exit_x(i)-person_x(i),exit_y(i)-person_y(i)].^2))*v0(i); %����i�������ٶ�����
                                Pl = Pl + (Kin*Vi_abs^h*abs(sum((Vj-Vi).*Ui))-C_ot)/max(0.2,Dij_abs-Radius(i)-Radius(j));
                            else %�н���0~30֮�䣬ͬʱ����Pl��Pm����
                                Vj = [vx(j),vy(j)]; %����j���ٶ�����Vj
                                Ui = [exit_x(i)-person_x(i),exit_y(i)-person_y(i)]/sqrt(sum([exit_x(i)-person_x(i),exit_y(i)-person_y(i)].^2))*v0(i); %����i�������ٶ�����
                                Pl = Pl + (Kin*Vi_abs^h*abs(sum((Vj-Vi).*Ui))-C_ot)/max(0.2,Dij_abs-Radius(i)-Radius(j));
                                Pm = Pm + (Kin*Vi_abs^h*abs(sum((Vj-Vi).*Ui))-C_ot)/max(0.2,Dij_abs-Radius(i)-Radius(j));
                            end
                        end
                    case 0 %VxD<0˵��Dij��Vi��˳ʱ�뷽�򣬼�j��i���Ҳ�
                        cos_VD = sum(Dij.*Vi)/(Dij_abs*Vi_abs);
                        if cos_VD<0.5 %�нǴ���60�㲻����
                            continue
                        else
                            if cos_VD<(sqrt(3)/2) %�н���30~60֮�䣬����Pr����
                                Vj = [vx(j),vy(j)]; %����j���ٶ�����Vj
                                Ui = [exit_x(i)-person_x(i),exit_y(i)-person_y(i)]/sqrt(sum([exit_x(i)-person_x(i),exit_y(i)-person_y(i)].^2))*v0(i); %����i�������ٶ�����
                                Pr = Pr + (Kin*Vi_abs^h*abs(sum((Vj-Vi).*Ui))-C_ot)/max(0.2,Dij_abs-Radius(i)-Radius(j));
                            else %�н���0~30֮�䣬ͬʱ����Pr��Pm����
                                Vj = [vx(j),vy(j)]; %����j���ٶ�����Vj
                                Ui = [exit_x(i)-person_x(i),exit_y(i)-person_y(i)]/sqrt(sum([exit_x(i)-person_x(i),exit_y(i)-person_y(i)].^2))*v0(i); %����i�������ٶ�����
                                Pr = Pr + (Kin*Vi_abs^h*abs(sum((Vj-Vi).*Ui))-C_ot)/max(0.2,Dij_abs-Radius(i)-Radius(j));
                                Pm = Pm + (Kin*Vi_abs^h*abs(sum((Vj-Vi).*Ui))-C_ot)/max(0.2,Dij_abs-Radius(i)-Radius(j));
                            end
                        end
                end
            end
        end
        S_l = w_p*Pl+w_rl*(C_rl-Vi_abs)*(-1);
        S_m = w_p*Pm+w_sa*d_sa*Vi_abs;
        S_r = w_p*Pr+w_rl*(C_rl-Vi_abs);
        index = find([S_l,S_m,S_r]==max([S_l,S_m,S_r]));
        Vi_0 = Vi/Vi_abs; %����i��ǰ�ٶȵĵ�λ����
        a = 60; %ת��Ƕ�
        switch index
            case 1 %���ֵΪS_l�������ĳ������ٶȷ���ΪVi��ʱ����ת90��
                a_pass_x(i) = a_pass_abs*(Vi_0(1)*cosd(a)-Vi_0(2)*sind(a));
                a_pass_y(i) = a_pass_abs*(Vi_0(1)*sind(a)+Vi_0(2)*cosd(a));
            case 2 %���ֵΪS_m���޳�����Ϊ
                continue
            case 3 %���ֵΪS_r�������ĳ������ٶȷ���ΪVi˳ʱ����ת90��
                a_pass_x(i) = a_pass_abs*(Vi_0(1)*cosd(-a)-Vi_0(2)*sind(-a));
                a_pass_y(i) = a_pass_abs*(Vi_0(1)*sind(-a)+Vi_0(2)*cosd(-a));
        end
    end
    
    %% �������˵�λ��
    ax = am_x+ar_x+ae_x+av_x+al_x+a_graX+a_pass_x;%1��n�У�tʱ�̸���������x����ĺϼ��ٶ�
    ay = am_y+ar_y+ae_y+av_y+al_y+a_graY+a_pass_y;%1��n�У�tʱ�̸���������y����ĺϼ��ٶ�
    %     ax = am_x+ar_x+ae_x+av_x+al_x+a_graX;%1��n�У�tʱ�̸���������x����ĺϼ��ٶ�
    %     ay = am_y+ar_y+ae_y+av_y+al_y+a_graY;%1��n�У�tʱ�̸���������y����ĺϼ��ٶ�
    vx = vx+ax*dt; %������һʱ�̵�x�����ٶ�
    vy = vy+ay*dt; %������һʱ�̵�y�����ٶ�
    V = sqrt(vx.^2+vy.^2);
    index = find(V>v0); %�ҳ��������ӵ�����
    vx(index) = vx(index).*v0(index)./V(index);
    vy(index) = vy(index).*v0(index)./V(index);
    
    person_x = person_x+vx*dt; %����x�����λ��
    person_y = person_y+vy*dt; %����y�����λ��
    
    %% ͳ�������ٶ�-�ܶ�ͼ
    %     index = find(person_x>49 & person_x<51); %Ѱ����x=46~52���䷶Χ�ڵ�����
    %     index_blue = find(index>person_l_num); %��x=46~52���䷶Χ�ڵ���ɫ����
    %     v_mean = mean(vx(index(index_blue))); %������ɫ�����ǵ�ƽ��ǰ���ٶ�
    %     if v_mean <= 0
    %         count = count+1; %�����ٶ�-�ܶ�ͳ��
    %         v_sum = v_sum+v_mean; %����������ƽ���ٶȺ�
    %         density_mean = length(index) / (2*4); %���������ܶ�
    %         density_sum = density_sum+density_mean; %���������������ܶ�֮��
    %     end
    %     if count == 20
    %         v_blue = [v_blue, v_sum/count];
    %         density_area = [density_area, density_sum/count];
    %         count = 0;
    %         v_sum = 0;
    %         density_sum = 0;
    %     end
    
    %% ͳ����ɢ��������Ĩ������ɢ���ӵ�������Ϣ
    for i=1:n
        if (person_x(i)-end_x(i))^2<=0.001
            sum_escape = sum_escape+1;
            person_x(i) = nan;
            person_y(i) = nan;
%             exit_x(i) = nan;
%             exit_y(i) = nan;
%             end_x(i) = nan;
%             vx(i) = nan;
%             vy(i) = nan;
%             v0(i) = nan;
%             Radius(i) = nan;
        end
    end
    %% ����ͼ��
    switch condition
        case 1
            % 4X50��ͼ
            index_l = find(exit_x==100);
            index_r = find(exit_x==-50);
            plot(wall_x1,wall_y1,'LineWidth',1,'Color','k');
            hold on;
            plot(wall_x2,wall_y2,'LineWidth',1,'Color','k');
            hold on;
            plot(person_x(index_l),person_y(index_l),'.r', 'MarkerSize', 10)
            plot(person_x(index_r),person_y(index_r),'.b', 'MarkerSize', 10)
            axis([0 50 -1 5]);%������ʾ��Χ
            if t==0 %ֻ�ڵ�1��ѭ������ͼ��λ�úʹ�С
                set(gcf,'position',[0,500,2000,260]);
            end
    end
    str_time=sprintf('��ɢʱ�䣺%.2f',t);
    str_escape=sprintf('��ɢ������%.0f',sum_escape);
    text(25,-0.5,str_time);
    text(25,-1.3,str_escape);
    axis on;
    hold off;
    pause(0.001);
    %% ɾ����Ĩ����Ϣ������
    index_del = find(isnan(person_x)==1);
    person_x(index_del) = [];
    person_y(index_del) = [];
    exit_x(index_del) = [];
    exit_y(index_del) = [];
    end_x(index_del) = [];
    vx(index_del) = [];
    vy(index_del) = [];
    v0(index_del) = [];
    Radius(index_del) = [];
end


