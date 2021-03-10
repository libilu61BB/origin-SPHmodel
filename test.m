% 尝试新的跟随模型（new-follow-model）

%% 更新日志
clear;
%% 设置障碍物坐标、行人坐标和出口坐标
condition = 1;
switch condition
    case 1
        % 4*100m通道及双向行人 17710813276 6-1-301
        personNum = 50; %每队行人的数量
        l_density = 1; %左侧行人的密度
        r_density = 1; %右侧行人的密度
        l_width = 40; %左侧行人初始化区域长度
        r_width = 40; %右侧行人初始化区域长度
        [person_l_num, person_l_x, person_l_y, person_r_num, person_r_x, person_r_y] = personInitialization(l_density,r_density,l_width,4,r_width,4);
        person_r_x = person_r_x + 60;
        person_x = [person_l_x person_r_x];
        person_y = [person_l_y person_r_y];
        wall_x1 = (0:0.1:100);
        wall_y1 = zeros(1, length(wall_x1));
        wall_x2 = (0:0.1:100);
        wall_y2 = 4 * ones(1, length(wall_x2));
        wall_x = [wall_x1, wall_x2];
        wall_y = [wall_y1, wall_y2];
%         person_x = [linspace(1,1+0.5*personNum,personNum), linspace(99-0.5*personNum,99,personNum)];
%         person_y = 3.4*rand(1,length(person_x))+0.3; %y∈[0.3 3.7]
        exit_x=[150*ones(1,person_l_num),-50*ones(1,person_r_num)];
        exit_y = 2*ones(1,person_l_num + person_r_num);
        end_x = [100*ones(1,person_l_num),0*ones(1,person_r_num)];
end
%% 计算坐标，绘制图像
n=length(person_x);
s=length(wall_x);
h1=5;%计算密度和排斥力时使用的核半径
Radius=0.3*ones(1,n);%假设行人的半径均为0.3m
m_person=70;%行人的质量
m_wall=500;%障碍物的质量
% va=1.2; %低速行人的期望速度值
% vb=1.8; %高速行人的期望速度值
v0 = [1.34*ones(1,floor(n/2)),1.34*ones(1,n-floor(n/2))]; %定义行人的期望速度，前一半为低速行人，后一半为高速行人,floor为向下取整
u=2; %粘度，用于计算粒子之间摩擦力产生的加速度
vx=zeros(1,n);%行人速度在x方向上的分量，初始时刻为0
vy=zeros(1,n);%行人速度在y方向上的分量，初始时刻为0
T=200; %模拟总时间
sum_escape=0;%统计已疏散的人数
P=1;%熟悉逃生路线的行人比例
P_f=1;%从众程度
dt=0.02;
% follow_A = 8;%跟随加速度的量级
tau = 0.2;%行人加速的特征时间

%% 朗之万随机力相关设置
P_r=0.5^dt;%加速度朗之万分量的时间权重
P_r = 0.5;%随便设一个试试
A=10;%加速度朗之万分量的量级
al=A*rand(1,n);%行人加速度的朗之万随机分量，为服从高斯分布的随机数
al_theta=2*pi*rand(1,n);%行人加速度朗之万随机分量的方向，为[0,2*pi]内的随机数
al_x=al.*cos(al_theta);%行人加速度在x方向上的朗之万随机力分量
al_y=al.*sin(al_theta);%行人加速度在y方向上的朗之万随机力分量

%% 速度-密度统计设置
count = 0; %速度密度统计计数
v_sum = 0; %平均速度和
density_sum = 0; %平均密度和
v_blue = [];
density_area = [];

%% 模拟循环
for t=0:dt:T
    %% 计算排斥力和挤压力产生的加速度
    [Rho_person,Rho_wall]=density(person_x,person_y,wall_x,wall_y,h1);%调用函数density计算t时刻的密度
    [ar_x,ar_y]=a_repul(person_x,person_y,wall_x,wall_y,Rho_person ,Rho_wall,Radius,h1);%调用函数a_repul计算排斥力产生的加速度
    [ae_x,ae_y]=a_extru(person_x,person_y,wall_x,wall_y,Rho_person,Rho_wall,Radius,h1);%调用函数a_extru计算挤压力产生的加速度
    %% 计算摩擦力产生的加速度
    av_x=zeros(1,n);%初始化摩擦力产生的x方向加速度
    av_y=zeros(1,n);%初始化摩擦力产生的y方向加速度
    for i=1:n
        for j=1:n %计算行人与行人之间的摩擦力产生的加速度
            if j==i
                continue;
            end
            r=sqrt((person_x(i)-person_x(j))^2+(person_y(i)-person_y(j))^2); %两行人粒子之间的距离
            if r<=(Radius(i)+Radius(j))
                av_x(i)=av_x(i)+u*m_person*((vx(j)-vx(i))/(Rho_person(i)*Rho_person(j)))*(1800*(Radius(i)+Radius(j)-r))/(45*pi*(Radius(i)+Radius(j))^5);
                av_y(i)=av_y(i)+u*m_person*((vy(j)-vy(i))/(Rho_person(i)*Rho_person(j)))*(1800*(Radius(i)+Radius(j)-r))/(45*pi*(Radius(i)+Radius(j))^5);
            end
        end
        d=zeros(1,s);%初始化行人粒子i与各障碍粒子的距离
        for j=1:s
            d(j)=sqrt((person_x(i)-wall_x(j))^2+(person_y(i)-wall_y(j))^2);
        end
        [r,m]=min(d); %r为d中最小值，m为d中最小值的索引
        if r<=Radius(i) %行人与障碍之间的摩擦力产生的加速度
            av_x(i)=av_x(i)+u*m_wall*((0-vx(i))/(Rho_person(i)*Rho_wall(m)))*(1800*(Radius(i)-r))/(45*pi*(Radius(i))^5);
            av_y(i)=av_y(i)+u*m_wall*((0-vy(i))/(Rho_person(i)*Rho_wall(m)))*(1800*(Radius(i)-r))/(45*pi*(Radius(i))^5);
        end
    end
    %% 计算行人主动力产生的加速度
    e_x=zeros(1,n); %初始化方向向量的x坐标
    e_y=zeros(1,n); %初始化方向向量的y坐标
    am_x=zeros(1,n); %初始化主动力产生的x方向加速度
    am_y=zeros(1,n); %初始化主动力产生的y方向加速度
    for i=1:(P*n) %计算熟悉逃生路线的行人的运动方向
        r=sqrt((person_x(i)-exit_x(i))^2+(person_y(i)-exit_y(i))^2); %熟悉逃生路线的行人粒子与出口之间的距离
        e_x(i)=(exit_x(i)-person_x(i))/r;
        e_y(i)=(exit_y(i)-person_y(i))/r;
    end
    for i=(P*n+1):n %为不熟悉逃生路线的行人产生一个随机的运动方向
        e_x(i)=-1+2*rand;
        e_y(i)=-1+2*rand;
        r=sqrt(e_x(i)^2+e_y(i)^2);%方向向量的模
        e_x(i)=e_x(i)/r;%将方向向量转换为单位方向向量
        e_y(i)=e_y(i)/r;%将方向向量转换为单位方向向量
    end
    for i=(P*n+1):n %计算不熟悉逃生路线的行人在从众行为影响下的运动方向
        e_x(i)=(1-P_f)*e_x(i);
        e_y(i)=(1-P_f)*e_y(i);
        for j=1:n
            if j==i
                continue;
            end
            r_2=(person_x(i)-person_x(j))^2+(person_y(i)-person_y(j))^2; %两行人粒子距离的平方
            if r_2<=h1^2
                e_x(i)=e_x(i)+P_f*(m_person*e_x(j)/Rho_person(j))*(4/(pi*h1^8))*(h1^2-r_2)^3;
                e_y(i)=e_y(i)+P_f*(m_person*e_y(j)/Rho_person(j))*(4/(pi*h1^8))*(h1^2-r_2)^3;
            end
        end
        r=sqrt(e_x(i)^2+e_y(i)^2);%方向向量的模
        e_x(i)=e_x(i)/r;%将方向向量转换为单位方向向量
        e_y(i)=e_y(i)/r;%将方向向量转换为单位方向向量
        r=sqrt((person_x(i)-exit_x(i))^2+(person_y(i)-exit_y(i))^2);
        if r<=5 %当不熟悉疏散路线的行人接近出口时，将其视为熟悉疏散路线的行人
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
    %% 计算朗之万随机加速度
    yita=A*rand(1,n);%行人加速度的朗之万随机分量，为服从高斯分布的随机数
    theta=2*pi*rand(1,n);%行人加速度朗之万随机分量的方向，为[0,2*pi]内的随机数
    al_x=P_r*al_x + (1-P_r)*yita.*cos(theta);%行人加速度在x方向上的朗之万随机力分量
    al_y=P_r*al_y + (1-P_r)*yita.*sin(theta);%行人加速度在y方向上的朗之万随机力分量
    %% 计算行人跟随行为产生的加速度
    % 根据引力模型的跟随行为
    a_graX = zeros(1,n); %初始化X跟随加速度
    a_graY = zeros(1,n); %初始化Y跟随加速度
    for i=1:n
        for j=1:n
            if i==j
                continue;
            end
            Dij = [person_x(j)-person_x(i),person_y(j)-person_y(i)]; %由i指向j的位置向量Dij
            Vi = [vx(i),vy(i)]; %粒子i的速度向量Vi
            Vj = [vx(j),vy(j)]; %粒子j的速度向量Vj
            L = sqrt(sum(Dij.^2)); %向量ij的模，相当于两粒子的距离
            if L<5 && sum(Dij.*Vi)>0 %当ij之间距离小于5且j位于i的前方时才有跟随行为，进行后续计算
                vj = sqrt(sum(Vj .^ 2)); %粒子j速度向量vj的模
                vi = sqrt(sum(Vi .^ 2)); %粒子i速度向量vi的模
                ej = Vj / vj; %粒子j速度方向向量ej
                ei0 = [exit_x(i)-person_x(i),exit_y(i)-person_y(i)]/sqrt(sum([exit_x(i)-person_x(i),exit_y(i)-person_y(i)].^2)); %由粒子i指向出口的单位向量
                a_graX(i) = a_graX(i) + v0(i)/tau * (vj-vi) * sum(ej .* ei0) / (L/(Radius(i)+Radius(j)))^2 * (person_x(j)-person_x(i))/L; %计算X跟随加速度
                a_graY(i) = a_graY(i) + v0(i)/tau * (vj-vi) * sum(ej .* ei0) / (L/(Radius(i)+Radius(j)))^2 * (person_y(j)-person_y(i))/L; %计算Y跟随加速度
            end
        end
    end
    
    
    
%     % 文献中的跟随行为模型
%     a_graX = zeros(1,n); %初始化X跟随加速度
%     a_graY = zeros(1,n); %初始化Y跟随加速度
%     for i=1:n
%         for j=1:n
%             if i==j
%                 continue;
%             end
%             Dij = [person_x(j)-person_x(i),person_y(j)-person_y(i)]; %由i指向j的位置向量Dij
%             Vi = [vx(i),vy(i)]; %粒子i的速度向量Vi
%             L = sqrt(sum(Dij.^2)); %向量ij的模，相当于两粒子的距离
%             if L<=2 && sum(Dij.*Vi)>0 %当距离小于等于2且j位于i的前方时才进行后续计算
%                 ei = [exit_x(i)-person_x(i),exit_y(i)-person_y(i)]/sqrt(sum([exit_x(i)-person_x(i),exit_y(i)-person_y(i)].^2)); %由粒子i指向出口的单位向量
%                 Vj = [vx(j),vy(j)]; %粒子j的速度向量Vj
%                 Bij_3 = max(0,sum(ei.*Vj)/sqrt(sum(Vj.^2)));
%                 if Bij_3==0
%                     continue;
%                 else
%                     Bij_4 = min(sqrt(sum(Vj.^2))/v0(i),1);
%                     Bij_5 = min(exp(1)^((Radius(i)+Radius(j))-L),1);
%                     a_graX(i) = a_graX(i)+0.4*v0(i)*Bij_3*Bij_4*Bij_5*Dij(1)/L; %计算X跟随加速度
%                     a_graY(i) = a_graY(i)+0.4*v0(i)*Bij_3*Bij_4*Bij_5*Dij(2)/L; %计算Y跟随加速度
%                 end         
%             end
%         end
%     end
    %% 计算行人的位置
%     ax = max(10,am_x+ar_x+ae_x+av_x+al_x+a_graX);%1行n列，t时刻各行人粒子x方向的合加速度
    ax = am_x+ar_x+ae_x+av_x+al_x+a_graX;%1行n列，t时刻各行人粒子x方向的合加速度
    ay = am_y+ar_y+ae_y+av_y+al_y+a_graY;%1行n列，t时刻各行人粒子y方向的合加速度
    vx = vx+ax*dt; %计算下一时刻的x方向速度
    vy = vy+ay*dt; %计算下一时刻的y方向速度
    V = sqrt(vx.^2+vy.^2);
    index = find(V>v0); %找出超速粒子的索引
    vx(index) = vx(index).*v0(index)./V(index);
    vy(index) = vy(index).*v0(index)./V(index);
    
    person_x = person_x+vx*dt; %计算x方向的位移
    person_y = person_y+vy*dt; %计算y方向的位移
    
    %% 统计行人速度-密度图
    
    index = find(person_x>49 & person_x<51); %寻找在x=46~52区间范围内的粒子
    index_blue = find(index>person_l_num); %在x=46~52区间范围内的蓝色粒子
    v_mean = mean(vx(index(index_blue))); %计算蓝色粒子们的平均前进速度
    if v_mean <= 0
        count = count+1; %进行速度-密度统计
        v_sum = v_sum+v_mean; %几个步长内平均速度和
        density_mean = length(index) / (2*4); %计算区域密度
        density_sum = density_sum+density_mean; %几个步长内区域密度之和
    end
    if count == 20
        v_blue = [v_blue, v_sum/count];
        density_area = [density_area, density_sum/count];
        count = 0;
        v_sum = 0;
        density_sum = 0;
    end

    
    %% 绘制图像
    for i=1:n %统计逃离人数
        if (person_x(i)-end_x(i))^2<=0.001
            sum_escape=sum_escape+1;
            person_x(i)=nan;
            person_y(i)=nan;
        end
    end
    switch condition
        case 1
            % 2X100画图
            plot(wall_x1,wall_y1,'LineWidth',1,'Color','k');
            hold on;
            plot(wall_x2,wall_y2,'LineWidth',1,'Color','k');
            hold on;
            plot(person_x(1:person_l_num),person_y(1:person_l_num),'.r', 'MarkerSize', 10)
            plot(person_x(person_l_num+1:end),person_y(person_l_num+1:end),'.b', 'MarkerSize', 10)
            axis([-1 101 -1 5]);%设置显示范围
            set(gcf,'position',[0,500,2000,160]);
    end
    str_time=sprintf('疏散时间：%.2f',t);
    str_escape=sprintf('逃离人数：%.0f',sum_escape);
    text(5.5,-0.5,str_time);
    text(5.5,-1.3,str_escape);
    axis on;
    hold off;
    pause(0.001);
        if sum_escape>=n
            break;
        end
end
        
        
        
    
                
                