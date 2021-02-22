function P_e = P_extru(Rho,Rho_0)
%P_extru 计算粒子相互接触挤压产生的压强
%   Rho 各粒子的密度
%   Rho_0 粒子的临界密度
%   Radius 各行人粒子的半径
n=length(Rho);
K=1000; %常数
P_e=zeros(1,n);
for i=1:n
    P_e(i)=K*(Rho(i)-Rho_0)/Rho_0;
    if P_e(i)<0
        P_e(i)=0;%将计算结果为负数的压强修正为0
    end
end
end