% 行人随机生成，直至数量达到设置的上限，且不删除已离开的行人的数据
clear;
% 2021-03-10
% 引力跟随行为模型V2.1：j粒子的影响投影到i当前速度与期望速度的速度差方向
% 双向行人流超车行为模拟，4*50m通道及连续双向行人
% 2021-03-21
% 更改速度-密度统计方式，修改右侧空隙密度和超车加速度大小

% 2021-04-15
% 空隙位置改为以i为中心的5个空隙，空隙方向用期望速度来确定，角度分别取[60 30 0 -30 -60]
% 空隙密度是两个范围的交集，第一个范围是以i为中心，search_R为半径的半圆，第二个范围是以空隙为中心，h_k为半径的圆
% 空隙密度算法更改后，容易出现0密度，导致超车加速度非常大，所以还是选择计算搜索范围内所有粒子的密度
% 计算空隙密度时，添加列一项|vj-vi|,考虑行人的速度差，速度差越大，密度贡献越大
% 我发现加速度过小的主要原因是密度的数值太大，所以添加了空隙密度归一化，将空隙密度缩放到0~Rho_kLim之间

clear;
condition = 1; %选择模拟场景
%% 初始化参数
switch condition
    case 1 %随机生成行人，两股行人流相向而行
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
        Sum_person = 100; %行人粒子的上限
        n = Sum_person;
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
        person_x = [20 30]; %行人的x坐标
        person_y = [2 2]; %行人的y坐标
        exit_x = [100 -50]; %出口的x坐标
        exit_y = [2 2]; %出口的y坐标
        end_x = [50 0]; %清除点的坐标
        vx = [0 0]; %行人速度在x方向上的分量
        vy = [0 0]; %行人速度在y方向上的分量
        v0 = [1.2 1.2]; %行人的期望速度
        n=length(person_x);
    case 3 %两粒子同向而行
        % 障碍物相关参数
        wall_x1 = (-20:0.1:70);
        wall_y1 = zeros(1, length(wall_x1));
        wall_x2 = (-20:0.1:70);
        wall_y2 = 4 * ones(1, length(wall_x2));
        wall_x = [wall_x1, wall_x2];
        wall_y = [wall_y1, wall_y2];
        s=length(wall_x);
        % 行人相关参数
        person_x = [0 10]; %行人的x坐标
        person_y = [2 2]; %行人的y坐标
        exit_x = [100 100]; %出口的x坐标
        exit_y = [2 2]; %出口的y坐标
%         end_x = [50 50]; %清除点的坐标
        vx = [0 0]; %行人速度在x方向上的分量
        vy = [0 0]; %行人速度在y方向上的分量
        v0 = [2.4 1.2]; %行人的期望速度
        n=length(person_x);
    case 4
        % 4*100m通道及单向行人
        wall_x1 = (0:0.1:100);
        wall_y1 = zeros(1, length(wall_x1));
        wall_x2 = (0:0.1:100);
        wall_y2 = 4 * ones(1, length(wall_x2));
        wall_x = [wall_x1, wall_x2];
        wall_y = [wall_y1, wall_y2];
        person_x = linspace(1,1+0.5*100,50);
        person_y = 3.4*rand(1,length(person_x))+0.3; %y∈[0.3 3.7]
        n=length(person_x);
        s=length(wall_x);
        exit_x=150*ones(1,n);
        exit_y = 2*ones(1,n);
        end_x = 100*ones(1,n);
        vx=zeros(1,n);%行人速度在x方向上的分量，初始时刻为0
        vy=zeros(1,n);%行人速度在y方向上的分量，初始时刻为0
        v0 = [2*ones(1,floor(n/2)),1*ones(1,n-floor(n/2))]; %定义行人的期望速度，前一半为低速行人，后一半为高速行人,floor为向下取整
end

h1 = 3; %计算密度和排斥力时使用的核半径
m_person=70; %行人的质量
m_wall=500; %障碍物的质量
Radius = 0.3; %行人的半径
u=20; %粘度，用于计算粒子之间摩擦力产生的加速度
T=200; %模拟总时间
tau = 0.2;%行人加速的特征时间
sum_escape=0; %统计已疏散的人数
dt=0.02; %时间步长
t_person = 0.5; %生成粒子的时间间隔

disP2P = zeros(n); %两行人之间的距离，初始化为n阶方阵，i行j列表示行人ij之间的距离，是主对角线为0的对称阵
disP2W = zeros(s,n); %行人与障碍物的距离
disP2E = zeros(1,n); %行人与出口的距离

% 超车行为2.0相关参数
search_R = 5; %计算区域得分时的搜索半径
L_i2k = 0.6; %粒子i与空隙的距离
h_k = 6; %计算行人对空隙密度影响时使用的核半径
h_wk = 1; %计算障碍对空隙密度影响时使用的核半径
K_pass = 1; %超车加速度修成系数
K_foll = 1; %跟随加速度修正系数
Rho_kLim = 10; %空隙密度的区间[0,Rho_kLim]

% 朗之万随机力相关参数
P_r=0.5^dt;%加速度朗之万分量的时间权重
A=0.5;%加速度朗之万分量的量级

% 速度-密度统计设置
count = 0; %速度密度统计计数
v_sum = 0; %平均速度和
density_sum = 0; %平均密度和
v_blue = [];
density_area = [];

%% 模拟循环
for t=0:dt:T
    %% 随机生成行人粒子
    if condition==1 && mod(t,t_person)==0 && n<=Sum_person %每隔固定的时间随机生成一次，粒子总数不超过设置值
        person_x_temp = [0-5*rand 50+5*rand]; %在左右两侧各生成一个粒子
        person_y_temp = [0.3+3.4*rand 0.3+3.4*rand];
        person_x = [person_x person_x_temp];
        person_y = [person_y person_y_temp];
        exit_x = [exit_x 100 -50];
        exit_y = [exit_y 2 2];
%         end_x = [end_x 35 15];
        vx = [vx 0 0];
        vy = [vy 0 0];
        temp1 = rand;
        temp2 = rand;
        v0 = [v0 1.2*(temp1>0.5)+1.8*(temp1<0.5) 1.2*(temp2>0.5)+1.8*(temp2<0.5)]; %高速行人和低速行人随机分配
        % v0 = [v0 1.36 1.36];
        n=length(person_x);
    end

    %% 计算各种的距离
    for i=1:n
        %↓↓两行人粒子之间的距离↓↓
        for j=(i+1):n
            disP2P(i,j) = sqrt((person_x(i) - person_x(j))^2 + (person_y(i) - person_y(j))^2); 
            disP2P(j,i) = disP2P(i,j); %ij和ji的距离一样，直接复制
        end
        %↓↓行人与障碍之间的距离↓↓
        for j=1:s
            disP2W(j,i) = sqrt((person_x(i)-wall_x(j))^2+(person_y(i)-wall_y(j))^2);
        end
        %↓↓行人与出口的距离↓↓
        disP2E(i) = sqrt((person_x(i)-exit_x(i))^2+(person_y(i)-exit_y(i))^2);
    end
    
    %% 郎之万随机力相关计算
    al=A*rand(1,n);%行人加速度的朗之万随机分量，为服从高斯分布的随机数
    al_theta=2*pi*rand(1,n);%行人加速度朗之万随机分量的方向，为[0,2*pi]内的随机数
    al_x=al.*cos(al_theta);%行人加速度在x方向上的朗之万随机力分量
    al_y=al.*sin(al_theta);%行人加速度在y方向上的朗之万随机力分量
    
    %% 计算排斥力和挤压力产生的加速度
    [Rho_person,Rho_wall]=density(n, m_person, m_wall, h1, disP2P, disP2W);%调用函数density计算t时刻的密度
    [ar_x,ar_y]=a_repul(person_x,person_y,wall_x,wall_y,Rho_person,Rho_wall,Radius,disP2P,disP2W,m_person,m_wall,h1);%调用函数a_repul计算排斥力产生的加速度
    [ae_x,ae_y]=a_extru(person_x,person_y,wall_x,wall_y,Rho_person,Rho_wall,Radius,disP2P,disP2W,m_person,m_wall);%调用函数a_extru计算挤压力产生的加速度
    
    %% 计算摩擦力产生的加速度
    av_x = zeros(1,n);%初始化摩擦力产生的x方向加速度
    av_y = zeros(1,n);%初始化摩擦力产生的y方向加速度
    h_f = 2*Radius;
    for i=1:n
        % ↓↓行人之间的摩擦↓↓
        disp2p = disP2P(i,1:n); %粒子i与其他粒子的距离
        disp2p(i) = nan; %不考虑粒子与自己的距离
        index = find(disp2p<h_f); %找到发生摩擦的粒子的索引
        av_x(i) = sum(u*m_person*((vx(index)-vx(i))./(Rho_person(i)*Rho_person(index))).* ...
            (1800*(2*Radius-disp2p(index)))./(45*pi*(2*Radius)^5));
        av_y(i) = sum(u*m_person*((vy(index)-vy(i))./(Rho_person(i)*Rho_person(index))).* ...
            (1800*(2*Radius-disp2p(index)))/(45*pi*(2*Radius)^5));
    end
    % ↓↓行人与障碍物之间的摩擦↓↓
    [disp2w,ind_w] = min(disP2W(:,1:n)); %disP2W按列计算最小值，并返回每列最小值的行号，相当于障碍物粒子的索引
    ind_p = find(disp2w<Radius); %计算与障碍物发生摩擦的行人粒子的索引
    if ~isempty(ind_p) %如果有粒子与障碍物产生摩擦，则计算人与障碍物之间的摩擦加速度
        av_x(ind_p) = av_x(ind_p)+u*m_wall*((0-vx(ind_p))./(Rho_person(ind_p).*Rho_wall(ind_w(ind_p)))).*...
            (1800*(Radius-disp2w(ind_p)))/(45*pi*(Radius)^5);
        av_y(ind_p) = av_y(ind_p)+u*m_wall*((0-vy(ind_p))./(Rho_person(ind_p).*Rho_wall(ind_w(ind_p)))).*...
            (1800*(Radius-disp2w(ind_p)))/(45*pi*(Radius)^5);
    end

    %% 计算行人主动力产生的加速度
    e_x = (exit_x-person_x)./disP2E(1:n); %指向出口的单位向量，x向
    e_y = (exit_y-person_y)./disP2E(1:n); %指向出口的单位向量，y向
    ui0_x = e_x.*v0; %i的期望速度，x向
    ui0_y = e_y.*v0; %i的期望速度，y向
    am_x = (ui0_x-vx)/tau;
    am_y = (ui0_y-vy)/tau;

    %% 计算朗之万随机加速度
    yita=A*rand(1,n);%行人加速度的朗之万随机分量，为服从高斯分布的随机数
    theta=2*pi*rand(1,n);%行人加速度朗之万随机分量的方向，为[0,2*pi]内的随机数
    al_x=P_r*al_x + (1-P_r)*yita.*cos(theta);%行人加速度在x方向上的朗之万随机力分量
    al_y=P_r*al_y + (1-P_r)*yita.*sin(theta);%行人加速度在y方向上的朗之万随机力分量
    
    
    %% 超车和跟随，新
    a_graX = zeros(1,n);
    a_graY = zeros(1,n);
    a_pass_x = zeros(1,n);
    a_pass_y = zeros(1,n);
    for i=1:n
        if vx(i)==0 && vy(i)==0
            continue
        end
        disp2p = disP2P(i,1:n);
%         disp2p(i) = nan;
        ind_R1 = find(disp2p<search_R); %找到搜索半径内的粒子j
        ind_R1(ind_R1==i) = [];
        mark1 = zeros(1,n);
        mark2 = zeros(1,n);
        mark1(ind_R1) = (person_x(ind_R1)-person_x(i))*vx(i)+...
            (person_y(ind_R1)-person_y(i))*vy(i); %Dij与Vi的向量积，用于判断j在i的前方
        ind_R2 = find(mark1>0); %找到搜索半径内位于i前方的粒子j
        mark2(ind_R2) = (vx(ind_R2)-vx(i))*(ui0_x(i)-vx(i))+...
            (vy(ind_R2)-vy(i))*(ui0_y(i)-vy(i)); %(Vj-Vi)与(ui0-Vi)的内积判断
        ind_pass = find(mark2<0); %i想要超车的j
        ind_foll = find(mark2>0); %i想要跟随的j
        % 计算跟随加速度
        if ~isempty(ind_foll) %如果有要跟随的粒子，就计算跟随加速度
            eij_x = (person_x(ind_foll)-person_x(i))./disp2p(ind_foll);
            eij_y = (person_y(ind_foll)-person_y(i))./disp2p(ind_foll);
            a_graX(i) = K_foll*sum(mark2(ind_foll)/tau.*(eij_x*e_x(i)+eij_y*e_y(i))/...
                (disp2p(ind_foll)/(2*Radius)).^2.*eij_x);
            a_graY(i) = K_foll*sum(mark2(ind_foll)/tau.*(eij_x*e_x(i)+eij_y*e_y(i))/...
                (disp2p(ind_foll)/(2*Radius)).^2.*eij_y);
        end

        % 计算超车加速度
        if ~isempty(ind_pass) %如果有要超越的粒子，就计算超车加速度
            ind_search = [ind_pass ind_foll];
            alpha = [60 30 0 -30 -60]; %空隙与Vi的夹角           
%             vx_i = vx(i)*cosd(alpha)-vy(i)*sind(alpha); %把Vi逆时针旋转alpha度，用来确定5个空隙的方向
%             vy_i = vx(i)*sind(alpha)+vy(i)*cosd(alpha);
%             abs_vi = sqrt(vx_i.^2+vy_i.^2);
%             k_x = person_x(i)+L_i2k*vx_i./abs_vi; %空隙的坐标，x向，1*5
%             k_y = person_y(i)+L_i2k*vy_i./abs_vi; %空隙的坐标，y向，1*5           
            vx_k = e_x(i)*cosd(alpha)-e_y(i)*sind(alpha); %把Vi逆时针旋转alpha度，用来确定5个空隙的方向
            vy_k = e_x(i)*sind(alpha)+e_y(i)*cosd(alpha);
            k_x = person_x(i)+L_i2k*vx_k; %空隙的坐标，x向，1*5
            k_y = person_y(i)+L_i2k*vy_k; %空隙的坐标，y向，1*5           
            m = length(ind_search);
            tempk = ones(m,1); %临时矩阵，用于扩展矩阵k_x和k_y
            tempj = ones(1,5); %临时矩阵，用于扩展ind_pass对应的矩阵person_x和person_y
            k_xtemp = tempk*k_x; %矩阵乘法，将k_x按行复制成m行5列
            k_ytemp = tempk*k_y;
            person_xtemp = person_x(ind_search)'*tempj; %矩阵乘法，将person_x转置后按列复制成m行5列
            person_ytemp = person_y(ind_search)'*tempj;
            disK2J = sqrt((k_xtemp-person_xtemp).^2+(k_ytemp-person_ytemp).^2); %计算空隙与粒子j的距离，m行5列，每列储存一个空隙的距离信息
            RhoK2J = m_person*(4/(pi*h_k^8))*(h_k^2-disK2J.^2).^3; %计算j对空隙的密度贡献
            RhoK2J(RhoK2J<0) = 0;
            
            vx_j = vx(ind_search);
            vy_j = vy(ind_search);
            absVji = sqrt((vx_j-vx(i)).^2+(vy_j-vy(i)).^2); %|vj-vi|
            absVji = absVji'*ones(1,5); %扩展成m行5列的矩阵
            RhoK2J = RhoK2J.*absVji;
            
            if sum(size(RhoK2J))>=7 
                RhoK2J = sum(RhoK2J); %按列求和，得到1行5列的密度矩阵
            end
            disp2w = disP2W(:,i)'; %第1个行人与障碍粒子的距离，转置为行向量
            [~,ind_i2w] = min(disp2w); %找到与i最近的墙壁，返回索引
            disK2W = sqrt((k_x-wall_x(ind_i2w(1))).^2+(k_y-wall_y(ind_i2w(1))).^2); %计算空隙与障碍物的距离
            RhoK2W = m_wall*(4/(pi*h_wk^8))*(h_wk^2-disK2W.^2).^3; %计算障碍物对空隙的密度贡献
            RhoK2W(RhoK2W<0) = 0;
            Rho_K = RhoK2J+RhoK2W; %计算空隙的总密度
            
            Rho_K = Rho_kLim*Rho_K/sum(Rho_K);
            [~,ind_k] = min(Rho_K);
            K1 = sqrt(sum([ui0_x(i)-vx(i),ui0_y(i)-vy(i)].^2))/tau/L_i2k^2;
            % K2 = (vx_k*e_x(i)+vy_k*e_y(i));
            K2 = 1;
            a_pass_xtemp = K1*K2./(Rho_K(ind_k)+eps)*vx_k(ind_k);
            a_pass_ytemp = K1*K2./(Rho_K(ind_k)+eps)*vy_k(ind_k);
            
            a_pass_x(i) = K_pass*a_pass_xtemp;
            a_pass_y(i) = K_pass*a_pass_ytemp;
            
%             K1 = sqrt(sum([ui0_x(i)-vx(i),ui0_y(i)-vy(i)].^2))/tau/L_i2k^2;
% %             K2 = (vx_k*e_x(i)+vy_k*e_y(i));
%             K2 = 1;
%             a_pass_xtemp = K1*K2./(Rho_K+eps).*vx_k;
%             a_pass_ytemp = K1*K2./(Rho_K+eps).*vy_k;
%             
%             a_pass_xtemp(a_pass_xtemp>10) = 10;
%             a_pass_xtemp(a_pass_xtemp<-10) = -10;
%             a_pass_ytemp(a_pass_ytemp>10) = 10;
%             a_pass_ytemp(a_pass_ytemp<-10) = -10;
%             
%             a_pass_abs = sqrt(a_pass_xtemp.^2+a_pass_ytemp.^2);
%             a_pass_abs = [1.0 1.0 1.0 1.2 1.4].*a_pass_abs;
%             [~,ind_k] = max(a_pass_abs);
%             a_pass_x(i) = K_pass*a_pass_xtemp(ind_k(1));
%             a_pass_y(i) = K_pass*a_pass_ytemp(ind_k(1));
        end
    end
    
    %% 计算行人的位置
    % 超车加速度和跟随加速度分解到xy轴上
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
%     temp_index = find(person_x>22 & person_x<28 & exit_x==100); %寻找在区间范围内且期望往右运动的粒子
%     if length(temp_index)>2
%         if count == 0
%             index_area = find(person_x>22 & person_x<28); %寻找x在22~28区间范围内的粒子
%             index_right = find(person_x>22 & person_x<28 & exit_x==100); %寻找在区间范围内且期望往右运动的粒子
%             right_x0 = person_x(index_right); %记录统计初始时刻向右运动粒子坐标
%             count = count + 1;
%         elseif count == 20
%             right_xt = person_x(index_right); %记录统计的末时刻向右运动粒子坐标
%             v_mean = (right_xt-right_x0)/(count*dt); %计算统计的粒子在几个时间步长内的平均速度
%             if mean(v_mean) > 2 || mean(v_mean) <0.3
%                 a = mean(v_mean);
%             end
%             density_mean = length(index_area) / (6*4); %计算区域密度
%             count = 0;
%             v_blue = [v_blue, mean(v_mean)];
%             density_area = [density_area, density_mean];
%         else
%             count = count + 1;
%         end
%     end
    
    % ↓↓统计疏散人数，并抹掉已疏散粒子的所有信息↓↓
%     for i=1:n
%         if (person_x(i)-end_x(i))^2<=0.1
%             sum_escape = sum_escape+1;
%             person_x(i) = nan;
%             person_y(i) = nan;
%         end
%     end
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
        case 3 % 4X50画图，两粒子同向而行
            plot(wall_x1,wall_y1,'LineWidth',1,'Color','k');
            hold on;
            plot(wall_x2,wall_y2,'LineWidth',1,'Color','k');
            plot(person_x,person_y,'.r', 'MarkerSize', 10)
            plot(person_x,person_y,'.b', 'MarkerSize', 10)
            axis([0 50 -2 5]);%设置显示范围
            if t==0 %只在第1次循环设置图窗位置和大小
                set(gcf,'position',[0,500,2000,260]);
            end
        case 4
            % 4X100画图
            plot(wall_x1,wall_y1,'LineWidth',1,'Color','k');
            hold on;
            plot(wall_x2,wall_y2,'LineWidth',1,'Color','k');
            hold on;
            plot(person_x(1:25),person_y(1:25),'or', 'MarkerSize', 10)
            plot(person_x(26:50),person_y(26:50),'ob', 'MarkerSize', 10)
            axis([-1 101 -1 5]);%设置显示范围
            set(gcf,'position',[0,500,2000,160]);
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
%     %% 删除被抹掉信息的粒子
%     index_del = find(isnan(person_x)==1);
%     if ~isempty(index_del)
%         person_x(index_del) = [];
%         person_y(index_del) = [];
%         exit_x(index_del) = [];
%         exit_y(index_del) = [];
%         end_x(index_del) = [];
%         vx(index_del) = [];
%         vy(index_del) = [];
%         v0(index_del) = [];
%         Radius(index_del) = [];
%     end
end


