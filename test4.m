%% 更新日志
clear;
% 2021-03-10
% 引力跟随行为模型V2.1：j粒子的影响投影到i当前速度与期望速度的速度差方向
% 双向行人流超车行为模拟，4*50m通道及连续双向行人
clear;
condition = 1; %选择模拟场景
%% 初始化参数
wall_x1 = (-20:0.1:70);
wall_y1 = zeros(1, length(wall_x1));
wall_x2 = (-20:0.1:70);
wall_y2 = 4 * ones(1, length(wall_x2));
wall_x = [wall_x1, wall_x2];
wall_y = [wall_y1, wall_y2];
s=length(wall_x);

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
t_person = 0.1; %生成粒子的时间间隔

% 超车行为相关参数设置(参考文献)
% w_p = 0.3; %区域得分权重
% w_sa = 0.4; %直线前进的权重
% w_rl = 0.3; %左右超车或避让的权重
% search_R = 5; %计算区域得分时的搜索半径
% a_pass_abs = 30; %超车行为产生的加速度的大小
% d_sa = 3;
% C_rl = 1.5;
% C_ot = 1.25;
% Kin = -1;
% h = 4;
% a = 30; %转弯角度

% 超车行为2.0相关参数
search_R = 5; %计算区域得分时的搜索半径
L_j2k = 1; %粒子j与空隙的距离
h_k = 2; %计算空隙密度时使用的核半径

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

%% 朗之万随机力相关设置
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
    if mod(t,t_person)==0 %每隔固定的时间随机生成一次
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
                av_x(i) = av_x(i)+u*m_person*((vx(j)-vx(i))/(Rho_person(i)*Rho_person(j)))*(1800*(Radius(i)+Radius(j)-r))/(45*pi*(Radius(i)+Radius(j))^5);
                av_y(i) = av_y(i)+u*m_person*((vy(j)-vy(i))/(Rho_person(i)*Rho_person(j)))*(1800*(Radius(i)+Radius(j)-r))/(45*pi*(Radius(i)+Radius(j))^5);
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
    %% 计算行人跟随行为产生的加速度
    % 引力模型跟随行为V2.1——形式1：j粒子的影响投影到i当前速度与期望速度的速度差方向
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
                Vij = Vj - Vi; %粒子i与粒子j的速度向量差
                ei0 = [exit_x(i)-person_x(i),exit_y(i)-person_y(i)]/sqrt(sum([exit_x(i)-person_x(i),exit_y(i)-person_y(i)].^2)); %由粒子i指向出口的单位向量
                ui0 = v0(i) * ei0; %粒子i的期望速度向量
                eij = Dij/L;%由粒子i指向粒子j的方向向量
                a_graX(i) = a_graX(i) + sum(Vij.*(ui0-Vi))/tau * sum(eij .* ei0) / (L/(Radius(i)+Radius(j)))^2 * (person_x(j)-person_x(i))/L; %计算X跟随加速度
                a_graY(i) = a_graY(i) + sum(Vij.*(ui0-Vi))/tau * sum(eij .* ei0) / (L/(Radius(i)+Radius(j)))^2 * (person_y(j)-person_y(i))/L; %计算Y跟随加速度
            end
        end
    end
    
    %% 计算超车行为产生的加速度
    % 超车模型2.0
    a_pass_x = zeros(1,n);
    a_pass_y = zeros(1,n);
    for i=1:n
        k_x = [];
        k_y = [];
        Vi = [vx(i),vy(i)]; %粒子i的速度向量
        e_i0 = [exit_x(i)-person_x(i),exit_y(i)-person_y(i)]/sqrt(sum([exit_x(i)-person_x(i),exit_y(i)-person_y(i)].^2));
        Ui = e_i0*v0(i); %粒子i的期望速度向量
        % ↓↓计算空隙的坐标（空隙与粒子j的距离由L_j2k控制）↓↓
        for j=1:n
            if i==j
                continue
            end
            Dij = [person_x(j)-person_x(i), person_y(j)-person_y(i)]; %由i指向j的位置向量Dij
            Dij_abs = sqrt(sum(Dij.^2)); %ij之间的距离
            if Dij_abs<=search_R %只计算搜索半径内的粒子
                cosij = Dij(1)/Dij_abs;
                sinij = Dij(2)/Dij_abs;
                Dik_r = [Dij(1)+L_j2k*sinij , Dij(2)-L_j2k*cosij]; %ij右侧空隙的向量坐标
                Dik_l = [Dij(1)-L_j2k*sinij , Dij(2)+L_j2k*cosij]; %ij左侧空隙的向量坐标
                k_x = [k_x Dik_r(1) Dik_l(1)]; %所有空隙的x坐标
                k_y = [k_y Dik_r(2) Dik_l(2)]; %所有空隙的y坐标
            end          
        end
        if isempty(k_x) %如果没有空隙，则跳过后续计算
            continue
        end
        % ↓↓计算空隙密度↓↓
        n_ik = length(k_x); %需要计算的空隙数量
        Rho_k = zeros(1,n_ik); %初始化空隙密度
        for kk = 1:n_ik
            for j = 1:n
                if j==i %不考虑粒子i对空隙密度的影响
                    continue
                end
                r2=(k_x(kk)-person_x(j))^2+(k_y(kk)-person_y(j))^2;%计算粒子j与空隙kk之间的距离平方
                if r2<=h_k^2 %计算核半径范围内粒子j对空隙kk密度的影响
                    Rho_k(kk) = Rho_k(kk)+m_person*(4/(pi*h_k^8))*(h_k^2-r2)^3;
                end
                [d2_wk,~] = min((k_x(kk)-wall_x).^2+(k_y(kk)-wall_y).^2);%计算空隙与障碍的最小距离（平方）及索引
                if d2_wk<=h_k^2 %只计算核半径范围内最近障碍粒子对粒子i密度的影响
                    Rho_k(kk)=Rho_k(kk)+m_wall*(4/(pi*h_k^8))*(h_k^2-d2_wk)^3;
                end
            end
        end
        % ↓↓选取密度最小的空隙作为超车方向↓↓
        [Rho_min,ind] = min(Rho_k); %返回最小密度及其索引
        Dik = [k_x(ind)-person_x(i),k_y(ind)-person_y(i)]; %由粒子i指向空隙k的向量
        r_ik = sqrt(sum(Dik.^2)); %粒子i与空隙k的距离
        e_ik = Dik/r_ik;
        a_pass = sqrt(sum((Ui-Vi).^2))*sum(e_ik.*e_i0)/(tau*Rho_min*r_ik^2)*e_ik; %粒子i的超车加速度
        a_pass_x(i) = a_pass(1);
        a_pass_y(i) = a_pass(2);
    end
    
    
    
    
    %% 参考文献中的超车模型
%     Pl = 0; %左侧区域的得分
%     Pm = 0; %中间区域的得分
%     Pr = 0; %右侧区域的得分   
%     a_pass_x = zeros(1,n);
%     a_pass_y = zeros(1,n);    
%     for i=1:n
%         Vi = [vx(i),vy(i)]; %粒子i的速度向量
%         Vi_abs = sqrt(sum(Vi.^2)); %粒子i的速度大小
%         if Vi_abs==0
%             continue
%         end
%         for j=1:n
%             if i==j
%                 continue
%             end
%             Dij = [person_x(j)-person_x(i), person_y(j)-person_y(i)]; %由i指向j的位置向量Dij
%             Dij_abs = sqrt(sum(Dij.^2)); %ij之间的距离
%             if Dij_abs<=search_R
%                 VxD = Vi(1)*Dij(2)-Vi(2)*Dij(1); %Vi与Dij的叉乘
%                 switch (VxD>=0)
%                     case 1 %VxD>=0说明Dij在Vi的逆时针方向，即j在i的左侧
%                         cos_VD = sum(Dij.*Vi)/(Dij_abs*Vi_abs);
%                         if cos_VD<0.5 %夹角大于60°不计算
%                             continue
%                         else
%                             if cos_VD<(sqrt(3)/2) %夹角在30~60之间，属于Pl区域
%                                 Vj = [vx(j),vy(j)]; %粒子j的速度向量Vj
%                                 Ui = [exit_x(i)-person_x(i),exit_y(i)-person_y(i)]/sqrt(sum([exit_x(i)-person_x(i),exit_y(i)-person_y(i)].^2))*v0(i); %粒子i的期望速度向量
%                                 Pl = Pl + (Kin*Vi_abs^h*abs(sum((Vj-Vi).*Ui))-C_ot)/max(0.2,Dij_abs-Radius(i)-Radius(j));
%                             else %夹角在0~30之间，同时属于Pl和Pm区域
%                                 Vj = [vx(j),vy(j)]; %粒子j的速度向量Vj
%                                 Ui = [exit_x(i)-person_x(i),exit_y(i)-person_y(i)]/sqrt(sum([exit_x(i)-person_x(i),exit_y(i)-person_y(i)].^2))*v0(i); %粒子i的期望速度向量
%                                 Pl = Pl + (Kin*Vi_abs^h*abs(sum((Vj-Vi).*Ui))-C_ot)/max(0.2,Dij_abs-Radius(i)-Radius(j));
%                                 Pm = Pm + (Kin*Vi_abs^h*abs(sum((Vj-Vi).*Ui))-C_ot)/max(0.2,Dij_abs-Radius(i)-Radius(j));
%                             end
%                         end
%                     case 0 %VxD<0说明Dij在Vi的顺时针方向，即j在i的右侧
%                         cos_VD = sum(Dij.*Vi)/(Dij_abs*Vi_abs);
%                         if cos_VD<0.5 %夹角大于60°不计算
%                             continue
%                         else
%                             if cos_VD<(sqrt(3)/2) %夹角在30~60之间，属于Pr区域
%                                 Vj = [vx(j),vy(j)]; %粒子j的速度向量Vj
%                                 Ui = [exit_x(i)-person_x(i),exit_y(i)-person_y(i)]/sqrt(sum([exit_x(i)-person_x(i),exit_y(i)-person_y(i)].^2))*v0(i); %粒子i的期望速度向量
%                                 Pr = Pr + (Kin*Vi_abs^h*abs(sum((Vj-Vi).*Ui))-C_ot)/max(0.2,Dij_abs-Radius(i)-Radius(j));
%                             else %夹角在0~30之间，同时属于Pr和Pm区域
%                                 Vj = [vx(j),vy(j)]; %粒子j的速度向量Vj
%                                 Ui = [exit_x(i)-person_x(i),exit_y(i)-person_y(i)]/sqrt(sum([exit_x(i)-person_x(i),exit_y(i)-person_y(i)].^2))*v0(i); %粒子i的期望速度向量
%                                 Pr = Pr + (Kin*Vi_abs^h*abs(sum((Vj-Vi).*Ui))-C_ot)/max(0.2,Dij_abs-Radius(i)-Radius(j));
%                                 Pm = Pm + (Kin*Vi_abs^h*abs(sum((Vj-Vi).*Ui))-C_ot)/max(0.2,Dij_abs-Radius(i)-Radius(j));
%                             end
%                         end
%                 end
%             end
%         end
%         S_l = w_p*Pl+w_rl*(C_rl-Vi_abs)*(-1);
%         S_m = w_p*Pm+w_sa*d_sa*Vi_abs;
%         S_r = w_p*Pr+w_rl*(C_rl-Vi_abs);
%         [~,index] = max([S_l,S_m,S_r]); %获取最大值的索引
%         Vi_0 = Vi/Vi_abs; %粒子i当前速度的单位向量       
%         switch index
%             case 1 %最大值为S_l，产生的超车加速度方向为Vi逆时针旋转a°
%                 a_pass_x(i) = a_pass_abs*(Vi_0(1)*cosd(a)-Vi_0(2)*sind(a));
%                 a_pass_y(i) = a_pass_abs*(Vi_0(1)*sind(a)+Vi_0(2)*cosd(a));
%             case 2 %最大值为S_m，无超车行为
%                 continue
%             case 3 %最大值为S_r，产生的超车加速度方向为Vi顺时针旋转a°
%                 a_pass_x(i) = a_pass_abs*(Vi_0(1)*cosd(-a)-Vi_0(2)*sind(-a));
%                 a_pass_y(i) = a_pass_abs*(Vi_0(1)*sind(-a)+Vi_0(2)*cosd(-a));
%         end
%     end
    
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
        case 1
            % 4X50画图
            index_l = find(exit_x==100);
            index_r = find(exit_x==-50);
            plot(wall_x1,wall_y1,'LineWidth',1,'Color','k');
            hold on;
            plot(wall_x2,wall_y2,'LineWidth',1,'Color','k');
            hold on;
            plot(person_x(index_l),person_y(index_l),'.r', 'MarkerSize', 10)
            plot(person_x(index_r),person_y(index_r),'.b', 'MarkerSize', 10)
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


