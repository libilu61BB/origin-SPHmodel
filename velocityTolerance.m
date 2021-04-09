%% ������־
% �о�������Ϊ�У���ǰ�������ٶȵ����Ͷȶ���Ⱥ����̶�����ɢЧ�ʵ�Ӱ��
% �޸ĸ�����Ϊ�볬����Ϊ���жϷ���
% ������Ϊģ�ͣ�����ģ��
% ������Ϊģ�ͣ�j�Ŀ�϶����ģ��
clear;
condition = 1; %ѡ��ģ�ⳡ��
%% ��ʼ������
switch condition
    case 1 %�����������
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
        Radius = []; %���˵İ뾶
end

h1=5; %�����ܶȺ��ų���ʱʹ�õĺ˰뾶
m_person=70; %���˵�����
m_wall=500; %�ϰ��������
u=4; %ճ�ȣ����ڼ�������֮��Ħ���������ļ��ٶ�
T=200; %ģ����ʱ��
tau = 0.2;%���˼��ٵ�����ʱ��
sum_escape=0; %ͳ������ɢ������
P=1; %��Ϥ����·�ߵ����˱���
P_f=1; %���ڳ̶�
dt=0.02; %ʱ�䲽��
t_person = 0.5; %�������ӵ�ʱ����

% ������Ϊ2.0��ز���
search_R = 5; %��������÷�ʱ�������뾶
L_j2k = 1; %����j���϶�ľ���
h_k = 5; %�������˶Կ�϶�ܶ�Ӱ��ʱʹ�õĺ˰뾶
h_wk = 0.3; %�����ϰ��Կ�϶�ܶ�Ӱ��ʱʹ�õĺ˰뾶
a = 1; %iǰ����϶���ܶ�����ָ��

% ��֮���������ز���
P_r=0.5^dt;%���ٶ���֮�������ʱ��Ȩ��
A=10;%���ٶ���֮�����������

% ���������̶�ͳ����ز���
trace2r_index = []; % �����ߵĴ�ͳ������
trace2r_x = zeros(200, 1000);
trace2r_y = zeros(200, 1000);
trace2l_index = []; % �����ߵĴ�ͳ������
trace2l_x = zeros(200, 1000);
trace2l_y = zeros(200, 1000);



%% ģ��ѭ��
for t=0:dt:T
    %% ���������������
    if condition==1 && mod(t,t_person)==0 %ÿ���̶���ʱ���������һ��
        person_x_temp = [10-5*rand 40+5*rand]; %���������������һ������
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
%         v0 = [v0 1.2*(temp1>0.5)+1.8*(temp1<0.5) 1.2*(temp2>0.5)+1.8*(temp2<0.5)]; %�������˺͵��������������
        v0 = [v0 1.36 1.36];
        Radius = [Radius 0.3 0.3];
    end
%     if condition==4 && mod(t,t_person)==0 %ÿ���̶���ʱ���������һ��
%         person_x_temp = [-2*rand-2 -2*rand-2]; %���������������һ������
%         person_y_temp = [0.3+2.4*rand 0.3+2.4*rand];
%         person_x = [person_x person_x_temp];
%         person_y = [person_y person_y_temp];
%         exit_x = [exit_x 100 100];
%         exit_y = [exit_y 1.5 1.5];
%         end_x = [end_x 12 12];
%         vx = [vx 0 0];
%         vy = [vy 0 0];
%         temp1 = rand;
%         temp2 = rand;
% %         v0 = [v0 1.2*(temp1>0.5)+1.8*(temp1<0.5) 1.2*(temp2>0.5)+1.8*(temp2<0.5)]; %�������˺͵��������������
%         v0 = [v0 1.55+0.36*rand-0.36 1.55+0.36*rand-0.36];
%         Radius = [Radius 0.3 0.3];
%     end
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
                av_x(i) = av_x(i)+u*m_person*((vx(j)-vx(i))/...%̫���ˣ���һ��
                    (Rho_person(i)*Rho_person(j)))*(1800*(Radius(i)+Radius(j)-r))/(45*pi*(Radius(i)+Radius(j))^5);
                av_y(i) = av_y(i)+u*m_person*((vy(j)-vy(i))/...
                    (Rho_person(i)*Rho_person(j)))*(1800*(Radius(i)+Radius(j)-r))/(45*pi*(Radius(i)+Radius(j))^5);
            end
        end
        d=zeros(1,s);%��ʼ����������i����ϰ����ӵľ���
        for j=1:s
            d(j)=sqrt((person_x(i)-wall_x(j))^2+(person_y(i)-wall_y(j))^2);
        end
        [r,m]=min(d); %rΪd����Сֵ��mΪd����Сֵ������
        if r<=Radius(i) %�������ϰ�֮���Ħ���������ļ��ٶ�
            av_x(i)=av_x(i)+u*m_wall*((0-vx(i))/...
                (Rho_person(i)*Rho_wall(m)))*(1800*(Radius(i)-r))/(45*pi*(Radius(i))^5);
            av_y(i)=av_y(i)+u*m_wall*((0-vy(i))/...
                (Rho_person(i)*Rho_wall(m)))*(1800*(Radius(i)-r))/(45*pi*(Radius(i))^5);
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
            r2=(person_x(i)-person_x(j))^2+(person_y(i)-person_y(j))^2; %���������Ӿ����ƽ��
            if r2<=h1^2
                e_x(i)=e_x(i)+P_f*(m_person*e_x(j)/Rho_person(j))*(4/(pi*h1^8))*(h1^2-r2)^3;
                e_y(i)=e_y(i)+P_f*(m_person*e_y(j)/Rho_person(j))*(4/(pi*h1^8))*(h1^2-r2)^3;
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
    
    %% �������ͳ�����Ϊ�����ļ��ٶ�
    a_graX = zeros(1,n); %��ʼ��X������ٶ�
    a_graY = zeros(1,n); %��ʼ��Y������ٶ�
    a_pass_x = zeros(1,n);
    a_pass_y = zeros(1,n);
    for i=1:n
        k_x = [];
        k_y = [];
        Vi = [vx(i),vy(i)]; %����i���ٶ�����Vi
        Vi_abs = sqrt(sum(Vi.^2)); %Vi�Ĵ�С
        if Vi_abs==0 %��ֹ������û�и���ͳ�����Ϊ
            continue
        end
        ei0 = [exit_x(i)-person_x(i),exit_y(i)-person_y(i)]/...
            sqrt(sum([exit_x(i)-person_x(i),exit_y(i)-person_y(i)].^2)); %������iָ����ڵĵ�λ����
        ui0 = v0(i) * ei0; %����i�������ٶ�����
        index_follow = zeros(1,n); %i��Ҫ�����j
        index_pass = zeros(1,n); %i��Ҫ������j
        for j=1:n
            if i==j
                continue
            end
            Vj = [vx(j),vy(j)]; %����j���ٶ�����Vj
            Dij = [person_x(j)-person_x(i),person_y(j)-person_y(i)]; %��iָ��j��λ������Dij
            Dij_abs = sqrt(sum(Dij.^2)); %Dij�Ĵ�С
            flag1 = sum((Vi.*Dij))/Vi_abs/Dij_abs; %Vi��Dij�ļн�
%             flag2 = sum((Vj-Vi).*(ui0-Vi)); %�ڻ��ж�
            if Dij_abs<=search_R && flag1>0.5 %ֻ���������뾶���Ҽн�С��60�������
%                 if flag2>0 
                if sum(Vj.*(ui0-Vi)) > 1*sum(Vi.*(ui0-Vi))
                    index_follow(j) = j; %����Ҫ��������ӵ�����
                else
                    index_pass(j) = j; %����Ҫ��Խ�����ӵ�����
                end
            end          
        end
        index_follow(index_follow==0) = []; %ɾ������Ԫ��
        index_pass(index_pass==0) = []; %ɾ������Ԫ��
        n_follow = length(index_follow);
        n_pass = length(index_pass);
        %% ������ٶȼ���
        for k=1:n_follow
            j = index_follow(k); %ֻ������������������j
            Vj = [vx(j),vy(j)];
            Vij = Vj - Vi; %����i������j���ٶ�������
            Dij = [person_x(j)-person_x(i),person_y(j)-person_y(i)];
            Dij_abs = sqrt(sum(Dij.^2));
            eij = Dij/Dij_abs;
            a_graX(i) = a_graX(i) + 1*sum(Vij.*(ui0-Vi))/tau * sum(eij .* ei0) /...
                (Dij_abs/(Radius(i)+Radius(j)))^2 * eij(1); %����X������ٶ�
            a_graY(i) = a_graY(i) + 1*sum(Vij.*(ui0-Vi))/tau * sum(eij .* ei0) /...
                (Dij_abs/(Radius(i)+Radius(j)))^2 * eij(2); %����Y������ٶ�
        end
        %% �������ٶȼ��㣨��jΪ��׼���ɿ�϶��
        Dok_m = [person_x(i)+L_j2k*Vi(1)/Vi_abs,person_y(i)+L_j2k*Vi(2)/Vi_abs]; %iǰ����϶�����꣨������ԭ��oָ���϶k,��ͬ��
        Rho_k = zeros(1,2*n_pass+1); %��ʼ����϶�ܶ�
        % ��������iǰ���Ŀ�϶�ܶȡ�������Rho_k(1)�С�����
        for jj = 1:n
            r2_m = (Dok_m(1)-person_x(jj))^2+(Dok_m(2)-person_y(jj))^2;%��������j���϶m֮��ľ���ƽ��
            if r2_m<=h_k^2 %����˰뾶��Χ������j�Կ�϶kk�ܶȵ�Ӱ��
                Rho_k(1) = Rho_k(1)+a*m_person*(4/(pi*h_k^8))*(h_k^2-r2_m)^3;
            end
        end
        d2_m = min((Dok_m(1)-wall_x).^2+(Dok_m(2)-wall_y).^2);%�����϶���ϰ�����С���루ƽ����
        if d2_m<=h_wk^2 %ֻ����˰뾶��Χ������ϰ����Ӷ�����i�ܶȵ�Ӱ��
            Rho_k(1)=Rho_k(1)+a*m_wall*(4/(pi*h_wk^8))*(h_wk^2-d2_m)^3;
        end
        % ��������j���������϶���ܶȡ��Ҳ��϶�ܶȴ�����Rho_k(2*k)�У�����϶�ܶȴ�����Rho_k(2*k+1)�С�����
        for k=1:length(index_pass)
            j = index_pass(k); %ֻ������������������j
            Dij = [person_x(j)-person_x(i),person_y(j)-person_y(i)];
            Dij_abs = sqrt(sum(Dij.^2));
            cosij = Dij(1)/Dij_abs;
            sinij = Dij(2)/Dij_abs;
            Dok_r = [person_x(j)+L_j2k*sinij , person_y(j)-L_j2k*cosij]; %ij�Ҳ��϶������
            Dok_l = [person_x(j)-L_j2k*sinij , person_y(j)+L_j2k*cosij]; %ij����϶������
            % ���������������Ӷ�j���������϶�ܶȵ�Ӱ�����
            for jj = 1:n
                if jj==i %����������i�Կ�϶�ܶȵ�Ӱ��
                    continue
                end
                r2_r=(Dok_r(1)-person_x(jj))^2+(Dok_r(2)-person_y(jj))^2;%��������jj��j�Ҳ��϶r֮��ľ���ƽ��
                r2_l=(Dok_l(1)-person_x(jj))^2+(Dok_l(2)-person_y(jj))^2;%��������jj��j����϶l֮��ľ���ƽ��
                if r2_r<=h_k^2
                    Rho_k(2*k) = Rho_k(2*k)+m_person*(4/(pi*h_k^8))*(h_k^2-r2_r)^3 / 15;
                end
                if r2_l<=h_k^2
                    Rho_k(2*k+1) = Rho_k(2*k+1)+m_person*(4/(pi*h_k^8))*(h_k^2-r2_l)^3;
                end
            end
            % ���������ϰ����Ӷ�j���������϶�ܶȵ�Ӱ�����
            d2_r = min((Dok_r(1)-wall_x).^2+(Dok_r(2)-wall_y).^2);%����j�Ҳ��϶���ϰ�����С���루ƽ����
            d2_l = min((Dok_l(1)-wall_x).^2+(Dok_l(2)-wall_y).^2);%����j����϶���ϰ�����С���루ƽ����
            if d2_r<=h_k^2 
                Rho_k(2*k)=Rho_k(2*k)+m_wall*(4/(pi*h_k^8))*(h_k^2-d2_r)^3 / 15;
            end
            if d2_l<=h_wk^2
                Rho_k(2*k+1)=Rho_k(2*k+1)+m_wall*(4/(pi*h_wk^8))*(h_wk^2-d2_l)^3;
            end
        end
        % ��������ÿ����϶��i������������
        n_k = length(Rho_k);
        Dik = cell(1,n_k); %��iָ���϶��λ������������Ԫ�����飩
        eik = cell(1,n_k); %Dik�ĵ�λ����������Ԫ�����飩
        Dik_abs = zeros(1,n_k);%����ik�Ĵ�С
        a_pass_abs = zeros(1,n_k); %�������ٶȵĴ�С
        A_k = 100*sqrt(sum((ui0-Vi).^2)); %�������ٶ�ϵ��
        for k=1:n_k
            if k==1 %��϶��iǰ��
                Dok = Dok_m;
            else
                if mod(k,2)==0 %��϶��ĳһ��ij�������Ҳ�
                    ind_j = index_pass(k/2); %�ҵ��Ǹ�j��person_x�е�����
                    Dij = [person_x(ind_j)-person_x(i),person_y(ind_j)-person_y(i)];
                    Dij_abs = sqrt(sum(Dij.^2));
                    cosij = Dij(1)/Dij_abs;
                    sinij = Dij(2)/Dij_abs;
                    Dok = [person_x(ind_j)+L_j2k*sinij , person_y(ind_j)-L_j2k*cosij];
                else %��϶��ĳһ��ij���������
                    ind_j = index_pass(floor(k/2)); %�ҵ��Ǹ�j��person_x�е�����
                    Dij = [person_x(ind_j)-person_x(i),person_y(ind_j)-person_y(i)];
                    Dij_abs = sqrt(sum(Dij.^2));
                    cosij = Dij(1)/Dij_abs;
                    sinij = Dij(2)/Dij_abs;
                    Dok = [person_x(ind_j)-L_j2k*sinij , person_y(ind_j)+L_j2k*cosij];
                end
            end
            Dik{k} = [Dok(1)-person_x(i),Dok(2)-person_y(i)];
            Dik_abs(k) = sqrt(sum(Dik{k}.^2)); 
            eik{k} = Dik{k}/Dik_abs(k);
            a_pass_abs(k) = A_k*sum(eik{k}.*ei0)/(tau*Rho_k(k)*Dik_abs(k)^2); %����i�ĳ������ٶ�
        end
        [a_pass_max,ind] = max(a_pass_abs);
        a_pass = a_pass_max*eik{ind};
        a_pass_x(i) = a_pass(1);
        a_pass_y(i) = a_pass(2);
        %% �������ٶȼ��㣨��iΪ��׼���ɿ�϶��
        
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
    
    %% ͳ���ȶ����к����˹켣
    % ����ʱ��100s֮��ȡx=20~30�����ڵ�200�����ӹ켣����ͳ��
    if t > 50 && length(trace2r_index)<200 % ͳ�Ƶ���x=20���������ߵ�����
        temp_trace = find(person_x>19.5 & person_x<20.5 & exit_x==100); 
        for tp=1:length(temp_trace)
            flag_trace = ismember(temp_trace(tp), trace2r_index);
            if ~flag_trace
                trace2r_index = [trace2r_index temp_trace(tp)];
            end
        end
        trace2r_nan = isnan(trace2r_index);
        for tp=1:length(trace2r_index)
            if ~trace2r_nan(tp)
                if person_x(trace2r_index(tp))>=19.5 && person_x(trace2r_index(tp))<=30.5 % ���ӹ켣ͳ�Ƶ�x=30Ϊֹ
                    temp_array = trace2r_x(tp, :);
                    temp_count = sum(temp_array~=0);
                    trace2r_x(tp,temp_count+1) = person_x(trace2r_index(tp));
                    trace2r_y(tp,temp_count+1) = person_y(trace2r_index(tp));
                else
                    trace2r_index(tp)=nan;
                end
            end
        end
%         tp_del = find(isnan(trace2r_index)==1);
%         if ~isempty(tp_del)
%             trace2r_index(tp_del) = [];
%         end
    end
    if t > 50 && length(trace2l_index)<200 % ͳ�Ƶ���x=30���������ߵ�����
        temp_trace = find(person_x>29.5 & person_x<30.5 & exit_x==-50); 
        for tp=1:length(temp_trace)
            flag_trace = ismember(temp_trace(tp), trace2l_index);
            if ~flag_trace
                trace2l_index = [trace2l_index temp_trace(tp)];
            end
        end
        trace2l_nan = isnan(trace2l_index);
        for tp=1:length(trace2l_index)
            if ~trace2l_nan(tp)
                if person_x(trace2l_index(tp))>=19.5 && person_x(trace2l_index(tp))<=30.5 % ���ӹ켣ͳ�Ƶ�x=30Ϊֹ
                    temp_array = trace2l_x(tp, :);
                    temp_count = sum(temp_array~=0);
                    trace2l_x(tp,temp_count+1) = person_x(trace2l_index(tp));
                    trace2l_y(tp,temp_count+1) = person_y(trace2l_index(tp));
                else
                    trace2l_index(tp)=nan;
                end
            end
        end
%         tp_del = find(isnan(trace2l_index)==1);
%         if ~isempty(tp_del)
%             trace2l_index(tp_del) = [];
%         end
    end

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
    end
    str_time = sprintf('��ɢʱ�䣺%.2f',t);
    str_escape = sprintf('��ɢ������%.0f',sum_escape);
    str_person = sprintf('��ǰ��������%d',n);
    text(5,-0.5,str_time);
    text(5,-1,str_escape);
    text(5,-1.5,str_person);
    axis on;
    hold off;
    pause(0.001);
    %% ɾ����Ĩ����Ϣ������
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

xlswrite('2rx_tolerance3.xlsx', trace2r_x);
xlswrite('2ry_tolerance3.xlsx', trace2r_y);
xlswrite('2lx_tolerance3.xlsx', trace2l_x);
xlswrite('2ly_tolerance3.xlsx', trace2l_y);



