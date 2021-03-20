%% 更新日志
clear;
% 2021-03-10
% 引力跟随行为模型V2.1：j粒子的影响投影到i当前速度与期望速度的速度差方向
% 双向行人流超车行为模拟，4*50m通道及连续双向行人
clear;
condition = 1; %选择模拟场景
%% 初始化参数
switch condition
    case 1 %随机生成行人
        % 障碍物相关参数
        wall_x1 = (-20:0.1:70);
        wall_y1 = zeros(1, length(wall_x1));
        wall_x2 = (-20:0.1:70);
        wall_y2 = 4 * ones(1, length(wall_x2));
        wall_x = [wall_x1, wall_x2];
        wall_y = [wall_y1, wall_y2];
        s=length(wall_x);
        % 行人相关参数
        person_x = []; %行人的x坐标
        person_y = []; %行人的y坐标
        exit_x = []; %出口的x坐标
        exit_y = []; %出口的y坐标
        end_x = []; %清除点的坐标
        vx = []; %行人速度在x方向上的分量
        vy = []; %行人速度在y方向上的分量
        v0 = []; %行人的期望速度
        Radius = []; %行人的半径
    case 2 %两粒子相向而行
        % 障碍物相关参数
        wall_x1 = (-20:0.1:70);
        wall_y1 = zeros(1, length(wall_x1));
        wall_x2 = (-20:0.1:70);
        wall_y2 = 4 * ones(1, length(wall_x2));
        wall_x = [wall_x1, wall_x2];
        wall_y = [wall_y1, wall_y2];
        s=length(wall_x);
        % 行人相关参数
        person_x = [10 40]; %行人的x坐标
        person_y = [2 2]; %行人的y坐标
        exit_x = [100 -50]; %出口的x坐标
        exit_y = [2 2]; %出口的y坐标
        end_x = [50 0]; %清除点的坐标
        vx = [0 0]; %行人速度在x方向上的分量
        vy = [0 0]; %行人速度在y方向上的分量
        v0 = [1.2 1.2]; %行人的期望速度
        Radius = [0.3 0.3]; %行人的半径       
end

h1=5; %计算密度和排斥力时使用的核半径
m_person=70; %行人的质量
m_wall=500; %障碍物的质量
u=2; %粘度，用于计算粒子之间摩擦力产生的加速度
T=200; %模拟总时间
tau = 0.2;%行人加速的特征时间
sum_escape=0; %统计已疏散的人数
P=1; %熟悉逃生路线的行人比例
P_f=1; %从众程度
dt=0.02; %时间步长
t_person = 0.5; %生成粒子的时间间隔

% 超车行为2.0相关参数
search_R = 5; %计算区域得分时的搜索半径
L_j2k = 0.6; %粒子j与空隙的距离
h_k = 2; %计算空隙密度时使用的核半径
a = 0.8; %i前方空隙的密度修正指数

% 朗之万随机力相关参数
P_r=0.5^dt;%加速度朗之万分量的时间权重
A=10;%加速度朗之万分量的量级

%% 速度-密度统计设置
count = 0; %速度密度统计计数
v_sum = 0; %平均速度和
density_sum = 0; %平均密度和
v_blue = [];
density_area = [];

%% 模拟循环
for t=0:dt:T
    %% 随机生成行人粒子
    if condition==1 && mod(t,t_person)==0 %每隔固定的时间随机生成一次
        person_x_temp = [-0.5-5*rand 50.5+5*rand]; %在左右两侧各生成一个粒子
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
        v0 = [v0 1.2*(temp1>0.5)+1.8*(temp1<0.5) 1.2*(temp2>0.5)+1.8*(temp2<0.5)]; %高速行人和低速行人随机分配
        Radius = [Radius 0.3 0.3];
    end
    n=length(person_x);
    %% 郎之万随机力相关计算
    al=A*rand(1,n);%行人加速度的朗之万随机分量，为服从高斯分布的随机数
    al_theta=2*pi*rand(1,n);%行人加速度朗之万随机分量的方向，为[0,2*pi]内的随机数
    al_x=al.*cos(al_theta);%行人加速度在x方向上的朗之万随机力分量
    al_y=al.*sin(al_theta);%行人加速度在y方向上的朗之万随机力分量
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
            r=sqrt((person_x(i) - person_x(j))^2 + (person_y(i) - person_y(j))^2); %两行人粒子之间的距离
            if r<=(Radius(i)+Radius(j))
                av_x(i) = av_x(i)+u*m_person*((vx(j)-vx(i))/...%太长了，换一行
                    (Rho_person(i)*Rho_person(j)))*(1800*(Radius(i)+Radius(j)-r))/(45*pi*(Radius(i)+Radius(j))^5);
                av_y(i) = av_y(i)+u*m_person*((vy(j)-vy(i))/...
                    (Rho_person(i)*Rho_person(j)))*(1800*(Radius(i)+Radius(j)-r))/(45*pi*(Radius(i)+Radius(j))^5);
            end
        end
        d=zeros(1,s);%初始化行人粒子i与各障碍粒子的距离
        for j=1:s
            d(j)=sqrt((person_x(i)-wall_x(j))^2+(person_y(i)-wall_y(j))^2);
        end
        [r,m]=min(d); %r为d中最小值，m为d中最小值的索引
        if r<=Radius(i) %行人与障碍之间的摩擦力产生的加速度
            av_x(i)=av_x(i)+u*m_wall*((0-vx(i))/...
                (Rho_person(i)*Rho_wall(m)))*(1800*(Radius(i)-r))/(45*pi*(Radius(i))^5);
            av_y(i)=av_y(i)+u*m_wall*((0-vy(i))/...
                (Rho_person(i)*Rho_wall(m)))*(1800*(Radius(i)-r))/(45*pi*(Radius(i))^5);
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
            r2=(person_x(i)-person_x(j))^2+(person_y(i)-person_y(j))^2; %两行人粒子距离的平方
            if r2<=h1^2
                e_x(i)=e_x(i)+P_f*(m_person*e_x(j)/Rho_person(j))*(4/(pi*h1^8))*(h1^2-r2)^3;
                e_y(i)=e_y(i)+P_f*(m_person*e_y(j)/Rho_person(j))*(4/(pi*h1^8))*(h1^2-r2)^3;
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
    
    %% 计算跟随和超车行为产生的加速度
    a_graX = zeros(1,n); %初始化X跟随加速度
    a_graY = zeros(1,n); %初始化Y跟随加速度
    a_pass_x = zeros(1,n);
    a_pass_y = zeros(1,n);
    for i=1:n
        k_x = [];
        k_y = [];
        Vi = [vx(i),vy(i)]; %粒子i的速度向量Vi
        Vi_abs = sqrt(sum(Vi.^2)); %Vi的大小
        if Vi_abs==0 %静止的粒子没有跟随和超车行为
            continue
        end
        ei0 = [exit_x(i)-person_x(i),exit_y(i)-person_y(i)]/...
            sqrt(sum([exit_x(i)-person_x(i),exit_y(i)-person_y(i)].^2)); %由粒子i指向出口的单位向量
        ui0 = v0(i) * ei0; %粒子i的期望速度向量
        index_follow = zeros(1,n); %i想要跟随的j
        index_pass = zeros(1,n); %i想要超车的j
        for j=1:n
            if i==j
                continue
            end
            Vj = [vx(j),vy(j)]; %粒子j的速度向量Vj
            Dij = [person_x(j)-person_x(i),person_y(j)-person_y(i)]; %由i指向j的位置向量Dij
            Dij_abs = sqrt(sum(Dij.^2)); %Dij的大小
            flag1 = sum((Vi.*Dij))/Vi_abs/Dij_abs; %Vi与Dij的夹角
            flag2 = sum((Vj-Vi).*(ui0-Vi)); %内积判断
            if Dij_abs<=search_R && flag1>0.5 %只考虑搜索半径内且夹角小于60°的粒子
                if flag2>0 
                    index_follow(j) = j; %储存要跟随的粒子的索引
                else
                    index_pass(j) = j; %储存要超越的粒子的索引
                end
            end          
        end
        index_follow(index_follow==0) = []; %删除多余元素
        index_pass(index_pass==0) = []; %删除多余元素
        n_follow = length(index_follow);
        n_pass = length(index_pass);
        %% 跟随加速度计算
        for k=1:n_follow
            j = index_follow(k); %只计算满足条件的粒子j
            Vj = [vx(j),vy(j)];
            Vij = Vj - Vi; %粒子i与粒子j的速度向量差
            Dij = [person_x(j)-person_x(i),person_y(j)-person_y(i)];
            Dij_abs = sqrt(sum(Dij.^2));
            eij = Dij/Dij_abs;
            a_graX(i) = a_graX(i) + sum(Vij.*(ui0-Vi))/tau * sum(eij .* ei0) /...
                (Dij_abs/(Radius(i)+Radius(j)))^2 * eij(1); %计算X跟随加速度
            a_graY(i) = a_graY(i) + sum(Vij.*(ui0-Vi))/tau * sum(eij .* ei0) /...
                (Dij_abs/(Radius(i)+Radius(j)))^2 * eij(2); %计算Y跟随加速度
        end
        %% 超车加速度计算
        Dok_m = [person_x(i)+L_j2k*Vi(1)/Vi_abs,person_y(i)+L_j2k*Vi(2)/Vi_abs]; %i前方空隙的坐标（由坐标原点o指向空隙k,下同）
        Rho_k = zeros(1,2*n_pass+1); %初始化空隙密度
        % ↓↓计算i前方的空隙密度【储存在Rho_k(1)中】↓↓
        for jj = 1:n
            r2_m = (Dok_m(1)-person_x(jj))^2+(Dok_m(2)-person_y(jj))^2;%计算粒子j与空隙m之间的距离平方
            if r2_m<=h_k^2 %计算核半径范围内粒子j对空隙kk密度的影响
                Rho_k(1) = Rho_k(1)+a*m_person*(4/(pi*h_k^8))*(h_k^2-r2_m)^3;
            end
        end
        d2_m = min((Dok_m(1)-wall_x).^2+(Dok_m(2)-wall_y).^2);%计算空隙与障碍的最小距离（平方）
        if d2_m<=h_k^2 %只计算核半径范围内最近障碍粒子对粒子i密度的影响
            Rho_k(1)=Rho_k(1)+a*m_wall*(4/(pi*h_k^8))*(h_k^2-d2_m)^3;
        end
        % ↓↓计算j左右两侧空隙的密度【右侧空隙密度储存在Rho_k(2*k)中，左侧空隙密度储存在Rho_k(2*k+1)中】↓↓
        for k=1:length(index_pass)
            j = index_pass(k); %只计算满足条件的粒子j
            Dij = [person_x(j)-person_x(i),person_y(j)-person_y(i)];
            Dij_abs = sqrt(sum(Dij.^2));
            cosij = Dij(1)/Dij_abs;
            sinij = Dij(2)/Dij_abs;
            Dok_r = [person_x(j)+L_j2k*sinij , person_y(j)-L_j2k*cosij]; %ij右侧空隙的坐标
            Dok_l = [person_x(j)-L_j2k*sinij , person_y(j)+L_j2k*cosij]; %ij左侧空隙的坐标
            % ↓↓计算行人粒子对j左右两侧空隙密度的影响↓↓
            for jj = 1:n
                if jj==i %不考虑粒子i对空隙密度的影响
                    continue
                end
                r2_r=(Dok_r(1)-person_x(jj))^2+(Dok_r(2)-person_y(jj))^2;%计算粒子jj与j右侧空隙r之间的距离平方
                r2_l=(Dok_l(1)-person_x(jj))^2+(Dok_l(2)-person_y(jj))^2;%计算粒子jj与j左侧空隙l之间的距离平方
                if r2_r<=h_k^2
                    Rho_k(2*k) = Rho_k(2*k)+m_person*(4/(pi*h_k^8))*(h_k^2-r2_r)^3;
                end
                if r2_l<=h_k^2
                    Rho_k(2*k+1) = Rho_k(2*k+1)+m_person*(4/(pi*h_k^8))*(h_k^2-r2_l)^3;
                end
            end
            % ↓↓计算障碍粒子对j左右两侧空隙密度的影响↓↓
            d2_r = min((Dok_r(1)-wall_x).^2+(Dok_r(2)-wall_y).^2);%计算j右侧空隙与障碍的最小距离（平方）
            d2_l = min((Dok_l(1)-wall_x).^2+(Dok_l(2)-wall_y).^2);%计算j左侧空隙与障碍的最小距离（平方）
            if d2_r<=h_k^2 
                Rho_k(2*k)=Rho_k(2*k)+m_wall*(4/(pi*h_k^8))*(h_k^2-d2_r)^3;
            end
            if d2_l<=h_k^2
                Rho_k(2*k+1)=Rho_k(2*k+1)+m_wall*(4/(pi*h_k^8))*(h_k^2-d2_l)^3;
            end
        end    
        % ↓↓选取密度最小的空隙作为超车方向↓↓
        [Rho_min,ind] = min(Rho_k); %返回最小密度及其索引
        if ind==1 %空隙在i前方
            Dok = Dok_m;
        else
            if mod(ind,2)==0 %空隙在某一个ij向量的右侧
                ind_j = index_pass(floor(ind/2)); %找到那个j在person_x中的索引
                Dij = [person_x(ind_j)-person_x(i),person_y(ind_j)-person_y(i)];
                Dij_abs = sqrt(sum(Dij.^2));
                cosij = Dij(1)/Dij_abs;
                sinij = Dij(2)/Dij_abs;
                Dok = [person_x(ind_j)+L_j2k*sinij , person_y(ind_j)-L_j2k*cosij];
            else %空隙在某一个ij向量的左侧
                ind_j = index_pass(floor(ind/2)); %找到那个j在person_x中的索引
                Dij = [person_x(ind_j)-person_x(i),person_y(ind_j)-person_y(i)];
                Dij_abs = sqrt(sum(Dij.^2));
                cosij = Dij(1)/Dij_abs;
                sinij = Dij(2)/Dij_abs;
                Dok = [person_x(ind_j)-L_j2k*sinij , person_y(ind_j)+L_j2k*cosij];
            end
        end        
        Dik = [Dok(1)-person_x(i),Dok(2)-person_y(i)]; %由粒子i指向空隙k的向量
        Dik_abs = sqrt(sum(Dik.^2)); %粒子i与空隙k的距离
        e_ik = Dik/Dik_abs;
        a_pass = sqrt(sum((ui0-Vi).^2))*sum(e_ik.*ei0)/(tau*Rho_min*Dik_abs^2)*e_ik; %粒子i的超车加速度
        a_pass_x(i) = a_pass(1);
        a_pass_y(i) = a_pass(2);
    end 
    
    %% 计算行人的位置
    ax = am_x+ar_x+ae_x+av_x+al_x+a_graX+a_pass_x;%1行n列，t时刻各行人粒子x方向的合加速度
    ay = am_y+ar_y+ae_y+av_y+al_y+a_graY+a_pass_y;%1行n列，t时刻各行人粒子y方向的合加速度
    vx = vx+ax*dt; %计算下一时刻的x方向速度
    vy = vy+ay*dt; %计算下一时刻的y方向速度
    V = sqrt(vx.^2+vy.^2);
    index = find(V>v0); %找出超速粒子的索引
    vx(index) = vx(index).*v0(index)./V(index);
    vy(index) = vy(index).*v0(index)./V(index);
    person_x = person_x+vx*dt; %计算x方向的位移
    person_y = person_y+vy*dt; %计算y方向的位移
    
    %% 统计行人速度-密度图
    %     index = find(person_x>49 & person_x<51); %寻找在x=46~52区间范围内的粒子
    %     index_blue = find(index>person_l_num); %在x=46~52区间范围内的蓝色粒子
    %     v_mean = mean(vx(index(index_blue))); %计算蓝色粒子们的平均前进速度
    %     if v_mean <= 0
    %         count = count+1; %进行速度-密度统计
    %         v_sum = v_sum+v_mean; %几个步长内平均速度和
    %         density_mean = length(index) / (2*4); %计算区域密度
    %         density_sum = density_sum+density_mean; %几个步长内区域密度之和
    %     end
    %     if count == 20
    %         v_blue = [v_blue, v_sum/count];
    %         density_area = [density_area, density_sum/count];
    %         count = 0;
    %         v_sum = 0;
    %         density_sum = 0;
    %     end
    
    %% 统计疏散人数，并抹掉已疏散粒子的所有信息
    for i=1:n
        if (person_x(i)-end_x(i))^2<=0.001
            sum_escape = sum_escape+1;
            person_x(i) = nan;
            person_y(i) = nan;
        end
    end
    %% 绘制图像
    switch condition
        case 1 % 4X50画图，行人粒子随机生成          
            index_l = find(exit_x==100);
            index_r = find(exit_x==-50);
            plot(wall_x1,wall_y1,'LineWidth',1,'Color','k');
            hold on;
            plot(wall_x2,wall_y2,'LineWidth',1,'Color','k');
            plot(person_x(index_l),person_y(index_l),'.r', 'MarkerSize', 10)
            plot(person_x(index_r),person_y(index_r),'.b', 'MarkerSize', 10)
            axis([0 50 -2 5]);%设置显示范围
            if t==0 %只在第1次循环设置图窗位置和大小
                set(gcf,'position',[0,500,2000,260]);
            end
        case 2 % 4X50画图，两粒子相向而行
            plot(wall_x1,wall_y1,'LineWidth',1,'Color','k');
            hold on;
            plot(wall_x2,wall_y2,'LineWidth',1,'Color','k');
            plot(person_x,person_y,'.r', 'MarkerSize', 10)
            plot(person_x,person_y,'.b', 'MarkerSize', 10)
            axis([0 50 -2 5]);%设置显示范围
            if t==0 %只在第1次循环设置图窗位置和大小
                set(gcf,'position',[0,500,2000,260]);
            end           
    end
    str_time = sprintf('疏散时间：%.2f',t);
    str_escape = sprintf('疏散人数：%.0f',sum_escape);
    str_person = sprintf('当前粒子数：%d',n);
    text(25,-0.5,str_time);
    text(25,-1,str_escape);
    text(25,-1.5,str_person);
    axis on;
    hold off;
    pause(0.001);
    %% 删除被抹掉信息的粒子
    index_del = find(isnan(person_x)==1);
    if ~isempty(index_del)
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
end


