% ����������ɣ�ֱ�������ﵽ���õ����ޣ��Ҳ�ɾ�����뿪�����˵�����
clear;
% 2021-03-10
% ����������Ϊģ��V2.1��j���ӵ�Ӱ��ͶӰ��i��ǰ�ٶ��������ٶȵ��ٶȲ��
% ˫��������������Ϊģ�⣬4*50mͨ��������˫������
% 2021-03-21
% �����ٶ�-�ܶ�ͳ�Ʒ�ʽ���޸��Ҳ��϶�ܶȺͳ������ٶȴ�С

% 2021-04-15
% ��϶λ�ø�Ϊ��iΪ���ĵ�5����϶����϶�����������ٶ���ȷ�����Ƕȷֱ�ȡ[60 30 0 -30 -60]
% ��϶�ܶ���������Χ�Ľ�������һ����Χ����iΪ���ģ�search_RΪ�뾶�İ�Բ���ڶ�����Χ���Կ�϶Ϊ���ģ�h_kΪ�뾶��Բ
% ��϶�ܶ��㷨���ĺ����׳���0�ܶȣ����³������ٶȷǳ������Ի���ѡ�����������Χ���������ӵ��ܶ�
% �����϶�ܶ�ʱ�������һ��|vj-vi|,�������˵��ٶȲ�ٶȲ�Խ���ܶȹ���Խ��
% �ҷ��ּ��ٶȹ�С����Ҫԭ�����ܶȵ���ֵ̫����������˿�϶�ܶȹ�һ��������϶�ܶ����ŵ�0~Rho_kLim֮��

clear;
condition = 1; %ѡ��ģ�ⳡ��
%% ��ʼ������
switch condition
    case 1 %����������ˣ������������������
        % �ϰ�����ز���
        wall_x1 = (-20:0.1:70);
        wall_y1 = zeros(1, length(wall_x1));
        wall_x2 = (-20:0.1:70);
        wall_y2 = 4 * ones(1, length(wall_x2));
        wall_x = [wall_x1, wall_x2];
        wall_y = [wall_y1, wall_y2];
        s=length(wall_x);
        % ������ز���
        person_x = []; %���˵�x����
        person_y = []; %���˵�y����
        exit_x = []; %���ڵ�x����
        exit_y = []; %���ڵ�y����
        end_x = []; %����������
        vx = []; %�����ٶ���x�����ϵķ���
        vy = []; %�����ٶ���y�����ϵķ���
        v0 = []; %���˵������ٶ�
        Sum_person = 100; %�������ӵ�����
        n = Sum_person;
    case 2 %�������������
        % �ϰ�����ز���
        wall_x1 = (-20:0.1:70);
        wall_y1 = zeros(1, length(wall_x1));
        wall_x2 = (-20:0.1:70);
        wall_y2 = 4 * ones(1, length(wall_x2));
        wall_x = [wall_x1, wall_x2];
        wall_y = [wall_y1, wall_y2];
        s=length(wall_x);
        % ������ز���
        person_x = [20 30]; %���˵�x����
        person_y = [2 2]; %���˵�y����
        exit_x = [100 -50]; %���ڵ�x����
        exit_y = [2 2]; %���ڵ�y����
        end_x = [50 0]; %����������
        vx = [0 0]; %�����ٶ���x�����ϵķ���
        vy = [0 0]; %�����ٶ���y�����ϵķ���
        v0 = [1.2 1.2]; %���˵������ٶ�
        n=length(person_x);
    case 3 %������ͬ�����
        % �ϰ�����ز���
        wall_x1 = (-20:0.1:70);
        wall_y1 = zeros(1, length(wall_x1));
        wall_x2 = (-20:0.1:70);
        wall_y2 = 4 * ones(1, length(wall_x2));
        wall_x = [wall_x1, wall_x2];
        wall_y = [wall_y1, wall_y2];
        s=length(wall_x);
        % ������ز���
        person_x = [0 10]; %���˵�x����
        person_y = [2 2]; %���˵�y����
        exit_x = [100 100]; %���ڵ�x����
        exit_y = [2 2]; %���ڵ�y����
%         end_x = [50 50]; %����������
        vx = [0 0]; %�����ٶ���x�����ϵķ���
        vy = [0 0]; %�����ٶ���y�����ϵķ���
        v0 = [2.4 1.2]; %���˵������ٶ�
        n=length(person_x);
    case 4
        % 4*100mͨ������������
        wall_x1 = (0:0.1:100);
        wall_y1 = zeros(1, length(wall_x1));
        wall_x2 = (0:0.1:100);
        wall_y2 = 4 * ones(1, length(wall_x2));
        wall_x = [wall_x1, wall_x2];
        wall_y = [wall_y1, wall_y2];
        person_x = linspace(1,1+0.5*100,50);
        person_y = 3.4*rand(1,length(person_x))+0.3; %y��[0.3 3.7]
        n=length(person_x);
        s=length(wall_x);
        exit_x=150*ones(1,n);
        exit_y = 2*ones(1,n);
        end_x = 100*ones(1,n);
        vx=zeros(1,n);%�����ٶ���x�����ϵķ�������ʼʱ��Ϊ0
        vy=zeros(1,n);%�����ٶ���y�����ϵķ�������ʼʱ��Ϊ0
        v0 = [2*ones(1,floor(n/2)),1*ones(1,n-floor(n/2))]; %�������˵������ٶȣ�ǰһ��Ϊ�������ˣ���һ��Ϊ��������,floorΪ����ȡ��
end

h1 = 3; %�����ܶȺ��ų���ʱʹ�õĺ˰뾶
m_person=70; %���˵�����
m_wall=500; %�ϰ��������
Radius = 0.3; %���˵İ뾶
u=20; %ճ�ȣ����ڼ�������֮��Ħ���������ļ��ٶ�
T=200; %ģ����ʱ��
tau = 0.2;%���˼��ٵ�����ʱ��
sum_escape=0; %ͳ������ɢ������
dt=0.02; %ʱ�䲽��
t_person = 0.5; %�������ӵ�ʱ����

disP2P = zeros(n); %������֮��ľ��룬��ʼ��Ϊn�׷���i��j�б�ʾ����ij֮��ľ��룬�����Խ���Ϊ0�ĶԳ���
disP2W = zeros(s,n); %�������ϰ���ľ���
disP2E = zeros(1,n); %��������ڵľ���

% ������Ϊ2.0��ز���
search_R = 5; %��������÷�ʱ�������뾶
L_i2k = 0.6; %����i���϶�ľ���
h_k = 6; %�������˶Կ�϶�ܶ�Ӱ��ʱʹ�õĺ˰뾶
h_wk = 1; %�����ϰ��Կ�϶�ܶ�Ӱ��ʱʹ�õĺ˰뾶
K_pass = 1; %�������ٶ��޳�ϵ��
K_foll = 1; %������ٶ�����ϵ��
Rho_kLim = 10; %��϶�ܶȵ�����[0,Rho_kLim]

% ��֮���������ز���
P_r=0.5^dt;%���ٶ���֮�������ʱ��Ȩ��
A=0.5;%���ٶ���֮�����������

% �ٶ�-�ܶ�ͳ������
count = 0; %�ٶ��ܶ�ͳ�Ƽ���
v_sum = 0; %ƽ���ٶȺ�
density_sum = 0; %ƽ���ܶȺ�
v_blue = [];
density_area = [];

%% ģ��ѭ��
for t=0:dt:T
    %% ���������������
    if condition==1 && mod(t,t_person)==0 && n<=Sum_person %ÿ���̶���ʱ���������һ�Σ�������������������ֵ
        person_x_temp = [0-5*rand 50+5*rand]; %���������������һ������
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
        v0 = [v0 1.2*(temp1>0.5)+1.8*(temp1<0.5) 1.2*(temp2>0.5)+1.8*(temp2<0.5)]; %�������˺͵��������������
        % v0 = [v0 1.36 1.36];
        n=length(person_x);
    end

    %% ������ֵľ���
    for i=1:n
        %��������������֮��ľ������
        for j=(i+1):n
            disP2P(i,j) = sqrt((person_x(i) - person_x(j))^2 + (person_y(i) - person_y(j))^2); 
            disP2P(j,i) = disP2P(i,j); %ij��ji�ľ���һ����ֱ�Ӹ���
        end
        %�����������ϰ�֮��ľ������
        for j=1:s
            disP2W(j,i) = sqrt((person_x(i)-wall_x(j))^2+(person_y(i)-wall_y(j))^2);
        end
        %������������ڵľ������
        disP2E(i) = sqrt((person_x(i)-exit_x(i))^2+(person_y(i)-exit_y(i))^2);
    end
    
    %% ��֮���������ؼ���
    al=A*rand(1,n);%���˼��ٶȵ���֮�����������Ϊ���Ӹ�˹�ֲ��������
    al_theta=2*pi*rand(1,n);%���˼��ٶ���֮����������ķ���Ϊ[0,2*pi]�ڵ������
    al_x=al.*cos(al_theta);%���˼��ٶ���x�����ϵ���֮�����������
    al_y=al.*sin(al_theta);%���˼��ٶ���y�����ϵ���֮�����������
    
    %% �����ų����ͼ�ѹ�������ļ��ٶ�
    [Rho_person,Rho_wall]=density(n, m_person, m_wall, h1, disP2P, disP2W);%���ú���density����tʱ�̵��ܶ�
    [ar_x,ar_y]=a_repul(person_x,person_y,wall_x,wall_y,Rho_person,Rho_wall,Radius,disP2P,disP2W,m_person,m_wall,h1);%���ú���a_repul�����ų��������ļ��ٶ�
    [ae_x,ae_y]=a_extru(person_x,person_y,wall_x,wall_y,Rho_person,Rho_wall,Radius,disP2P,disP2W,m_person,m_wall);%���ú���a_extru���㼷ѹ�������ļ��ٶ�
    
    %% ����Ħ���������ļ��ٶ�
    av_x = zeros(1,n);%��ʼ��Ħ����������x������ٶ�
    av_y = zeros(1,n);%��ʼ��Ħ����������y������ٶ�
    h_f = 2*Radius;
    for i=1:n
        % ��������֮���Ħ������
        disp2p = disP2P(i,1:n); %����i���������ӵľ���
        disp2p(i) = nan; %�������������Լ��ľ���
        index = find(disp2p<h_f); %�ҵ�����Ħ�������ӵ�����
        av_x(i) = sum(u*m_person*((vx(index)-vx(i))./(Rho_person(i)*Rho_person(index))).* ...
            (1800*(2*Radius-disp2p(index)))./(45*pi*(2*Radius)^5));
        av_y(i) = sum(u*m_person*((vy(index)-vy(i))./(Rho_person(i)*Rho_person(index))).* ...
            (1800*(2*Radius-disp2p(index)))/(45*pi*(2*Radius)^5));
    end
    % �����������ϰ���֮���Ħ������
    [disp2w,ind_w] = min(disP2W(:,1:n)); %disP2W���м�����Сֵ��������ÿ����Сֵ���кţ��൱���ϰ������ӵ�����
    ind_p = find(disp2w<Radius); %�������ϰ��﷢��Ħ�����������ӵ�����
    if ~isempty(ind_p) %������������ϰ������Ħ��������������ϰ���֮���Ħ�����ٶ�
        av_x(ind_p) = av_x(ind_p)+u*m_wall*((0-vx(ind_p))./(Rho_person(ind_p).*Rho_wall(ind_w(ind_p)))).*...
            (1800*(Radius-disp2w(ind_p)))/(45*pi*(Radius)^5);
        av_y(ind_p) = av_y(ind_p)+u*m_wall*((0-vy(ind_p))./(Rho_person(ind_p).*Rho_wall(ind_w(ind_p)))).*...
            (1800*(Radius-disp2w(ind_p)))/(45*pi*(Radius)^5);
    end

    %% �������������������ļ��ٶ�
    e_x = (exit_x-person_x)./disP2E(1:n); %ָ����ڵĵ�λ������x��
    e_y = (exit_y-person_y)./disP2E(1:n); %ָ����ڵĵ�λ������y��
    ui0_x = e_x.*v0; %i�������ٶȣ�x��
    ui0_y = e_y.*v0; %i�������ٶȣ�y��
    am_x = (ui0_x-vx)/tau;
    am_y = (ui0_y-vy)/tau;

    %% ������֮��������ٶ�
    yita=A*rand(1,n);%���˼��ٶȵ���֮�����������Ϊ���Ӹ�˹�ֲ��������
    theta=2*pi*rand(1,n);%���˼��ٶ���֮����������ķ���Ϊ[0,2*pi]�ڵ������
    al_x=P_r*al_x + (1-P_r)*yita.*cos(theta);%���˼��ٶ���x�����ϵ���֮�����������
    al_y=P_r*al_y + (1-P_r)*yita.*sin(theta);%���˼��ٶ���y�����ϵ���֮�����������
    
    
    %% �����͸��棬��
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
        ind_R1 = find(disp2p<search_R); %�ҵ������뾶�ڵ�����j
        ind_R1(ind_R1==i) = [];
        mark1 = zeros(1,n);
        mark2 = zeros(1,n);
        mark1(ind_R1) = (person_x(ind_R1)-person_x(i))*vx(i)+...
            (person_y(ind_R1)-person_y(i))*vy(i); %Dij��Vi���������������ж�j��i��ǰ��
        ind_R2 = find(mark1>0); %�ҵ������뾶��λ��iǰ��������j
        mark2(ind_R2) = (vx(ind_R2)-vx(i))*(ui0_x(i)-vx(i))+...
            (vy(ind_R2)-vy(i))*(ui0_y(i)-vy(i)); %(Vj-Vi)��(ui0-Vi)���ڻ��ж�
        ind_pass = find(mark2<0); %i��Ҫ������j
        ind_foll = find(mark2>0); %i��Ҫ�����j
        % ���������ٶ�
        if ~isempty(ind_foll) %�����Ҫ��������ӣ��ͼ��������ٶ�
            eij_x = (person_x(ind_foll)-person_x(i))./disp2p(ind_foll);
            eij_y = (person_y(ind_foll)-person_y(i))./disp2p(ind_foll);
            a_graX(i) = K_foll*sum(mark2(ind_foll)/tau.*(eij_x*e_x(i)+eij_y*e_y(i))/...
                (disp2p(ind_foll)/(2*Radius)).^2.*eij_x);
            a_graY(i) = K_foll*sum(mark2(ind_foll)/tau.*(eij_x*e_x(i)+eij_y*e_y(i))/...
                (disp2p(ind_foll)/(2*Radius)).^2.*eij_y);
        end

        % ���㳬�����ٶ�
        if ~isempty(ind_pass) %�����Ҫ��Խ�����ӣ��ͼ��㳬�����ٶ�
            ind_search = [ind_pass ind_foll];
            alpha = [60 30 0 -30 -60]; %��϶��Vi�ļн�           
%             vx_i = vx(i)*cosd(alpha)-vy(i)*sind(alpha); %��Vi��ʱ����תalpha�ȣ�����ȷ��5����϶�ķ���
%             vy_i = vx(i)*sind(alpha)+vy(i)*cosd(alpha);
%             abs_vi = sqrt(vx_i.^2+vy_i.^2);
%             k_x = person_x(i)+L_i2k*vx_i./abs_vi; %��϶�����꣬x��1*5
%             k_y = person_y(i)+L_i2k*vy_i./abs_vi; %��϶�����꣬y��1*5           
            vx_k = e_x(i)*cosd(alpha)-e_y(i)*sind(alpha); %��Vi��ʱ����תalpha�ȣ�����ȷ��5����϶�ķ���
            vy_k = e_x(i)*sind(alpha)+e_y(i)*cosd(alpha);
            k_x = person_x(i)+L_i2k*vx_k; %��϶�����꣬x��1*5
            k_y = person_y(i)+L_i2k*vy_k; %��϶�����꣬y��1*5           
            m = length(ind_search);
            tempk = ones(m,1); %��ʱ����������չ����k_x��k_y
            tempj = ones(1,5); %��ʱ����������չind_pass��Ӧ�ľ���person_x��person_y
            k_xtemp = tempk*k_x; %����˷�����k_x���и��Ƴ�m��5��
            k_ytemp = tempk*k_y;
            person_xtemp = person_x(ind_search)'*tempj; %����˷�����person_xת�ú��и��Ƴ�m��5��
            person_ytemp = person_y(ind_search)'*tempj;
            disK2J = sqrt((k_xtemp-person_xtemp).^2+(k_ytemp-person_ytemp).^2); %�����϶������j�ľ��룬m��5�У�ÿ�д���һ����϶�ľ�����Ϣ
            RhoK2J = m_person*(4/(pi*h_k^8))*(h_k^2-disK2J.^2).^3; %����j�Կ�϶���ܶȹ���
            RhoK2J(RhoK2J<0) = 0;
            
            vx_j = vx(ind_search);
            vy_j = vy(ind_search);
            absVji = sqrt((vx_j-vx(i)).^2+(vy_j-vy(i)).^2); %|vj-vi|
            absVji = absVji'*ones(1,5); %��չ��m��5�еľ���
            RhoK2J = RhoK2J.*absVji;
            
            if sum(size(RhoK2J))>=7 
                RhoK2J = sum(RhoK2J); %������ͣ��õ�1��5�е��ܶȾ���
            end
            disp2w = disP2W(:,i)'; %��1���������ϰ����ӵľ��룬ת��Ϊ������
            [~,ind_i2w] = min(disp2w); %�ҵ���i�����ǽ�ڣ���������
            disK2W = sqrt((k_x-wall_x(ind_i2w(1))).^2+(k_y-wall_y(ind_i2w(1))).^2); %�����϶���ϰ���ľ���
            RhoK2W = m_wall*(4/(pi*h_wk^8))*(h_wk^2-disK2W.^2).^3; %�����ϰ���Կ�϶���ܶȹ���
            RhoK2W(RhoK2W<0) = 0;
            Rho_K = RhoK2J+RhoK2W; %�����϶�����ܶ�
            
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
    
    %% �������˵�λ��
    % �������ٶȺ͸�����ٶȷֽ⵽xy����
    ax = am_x+ar_x+ae_x+av_x+al_x+a_graX+a_pass_x;%1��n�У�tʱ�̸���������x����ĺϼ��ٶ�
    ay = am_y+ar_y+ae_y+av_y+al_y+a_graY+a_pass_y;%1��n�У�tʱ�̸���������y����ĺϼ��ٶ�  
    vx = vx+ax*dt; %������һʱ�̵�x�����ٶ�
    vy = vy+ay*dt; %������һʱ�̵�y�����ٶ�
    V = sqrt(vx.^2+vy.^2);
    index = find(V>v0); %�ҳ��������ӵ�����
    vx(index) = vx(index).*v0(index)./V(index);
    vy(index) = vy(index).*v0(index)./V(index);
    person_x = person_x+vx*dt; %����x�����λ��
    person_y = person_y+vy*dt; %����y�����λ��
    
    %% ͳ�������ٶ�-�ܶ�ͼ
%     temp_index = find(person_x>22 & person_x<28 & exit_x==100); %Ѱ�������䷶Χ�������������˶�������
%     if length(temp_index)>2
%         if count == 0
%             index_area = find(person_x>22 & person_x<28); %Ѱ��x��22~28���䷶Χ�ڵ�����
%             index_right = find(person_x>22 & person_x<28 & exit_x==100); %Ѱ�������䷶Χ�������������˶�������
%             right_x0 = person_x(index_right); %��¼ͳ�Ƴ�ʼʱ�������˶���������
%             count = count + 1;
%         elseif count == 20
%             right_xt = person_x(index_right); %��¼ͳ�Ƶ�ĩʱ�������˶���������
%             v_mean = (right_xt-right_x0)/(count*dt); %����ͳ�Ƶ������ڼ���ʱ�䲽���ڵ�ƽ���ٶ�
%             if mean(v_mean) > 2 || mean(v_mean) <0.3
%                 a = mean(v_mean);
%             end
%             density_mean = length(index_area) / (6*4); %���������ܶ�
%             count = 0;
%             v_blue = [v_blue, mean(v_mean)];
%             density_area = [density_area, density_mean];
%         else
%             count = count + 1;
%         end
%     end
    
    % ����ͳ����ɢ��������Ĩ������ɢ���ӵ�������Ϣ����
%     for i=1:n
%         if (person_x(i)-end_x(i))^2<=0.1
%             sum_escape = sum_escape+1;
%             person_x(i) = nan;
%             person_y(i) = nan;
%         end
%     end
    %% ����ͼ��
    switch condition
        case 1 % 4X50��ͼ�����������������          
            index_l = find(exit_x==100);
            index_r = find(exit_x==-50);
            plot(wall_x1,wall_y1,'LineWidth',1,'Color','k');
            hold on;
            plot(wall_x2,wall_y2,'LineWidth',1,'Color','k');
            plot(person_x(index_l),person_y(index_l),'.r', 'MarkerSize', 10)
            plot(person_x(index_r),person_y(index_r),'.b', 'MarkerSize', 10)
            axis([0 50 -2 5]);%������ʾ��Χ
            if t==0 %ֻ�ڵ�1��ѭ������ͼ��λ�úʹ�С
                set(gcf,'position',[0,500,2000,260]);
            end
        case 2 % 4X50��ͼ���������������
            plot(wall_x1,wall_y1,'LineWidth',1,'Color','k');
            hold on;
            plot(wall_x2,wall_y2,'LineWidth',1,'Color','k');
            plot(person_x,person_y,'.r', 'MarkerSize', 10)
            plot(person_x,person_y,'.b', 'MarkerSize', 10)
            axis([0 50 -2 5]);%������ʾ��Χ
            if t==0 %ֻ�ڵ�1��ѭ������ͼ��λ�úʹ�С
                set(gcf,'position',[0,500,2000,260]);
            end
        case 3 % 4X50��ͼ��������ͬ�����
            plot(wall_x1,wall_y1,'LineWidth',1,'Color','k');
            hold on;
            plot(wall_x2,wall_y2,'LineWidth',1,'Color','k');
            plot(person_x,person_y,'.r', 'MarkerSize', 10)
            plot(person_x,person_y,'.b', 'MarkerSize', 10)
            axis([0 50 -2 5]);%������ʾ��Χ
            if t==0 %ֻ�ڵ�1��ѭ������ͼ��λ�úʹ�С
                set(gcf,'position',[0,500,2000,260]);
            end
        case 4
            % 4X100��ͼ
            plot(wall_x1,wall_y1,'LineWidth',1,'Color','k');
            hold on;
            plot(wall_x2,wall_y2,'LineWidth',1,'Color','k');
            hold on;
            plot(person_x(1:25),person_y(1:25),'or', 'MarkerSize', 10)
            plot(person_x(26:50),person_y(26:50),'ob', 'MarkerSize', 10)
            axis([-1 101 -1 5]);%������ʾ��Χ
            set(gcf,'position',[0,500,2000,160]);
    end
    str_time = sprintf('��ɢʱ�䣺%.2f',t);
    str_escape = sprintf('��ɢ������%.0f',sum_escape);
    str_person = sprintf('��ǰ��������%d',n);
    text(25,-0.5,str_time);
    text(25,-1,str_escape);
    text(25,-1.5,str_person);
    axis on;
    hold off;
    pause(0.001);
%     %% ɾ����Ĩ����Ϣ������
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


