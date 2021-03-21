function P_r = P_repul(Rho,Rho_0)
%P_repul 计算排斥力带给每个行人的压强
%   Rho 各粒子的密度
%   Rho_0 粒子的临界密度
%   A 压强计算系数
%   B 压强计算系数
n=length(Rho);
A=125;
B=0.001;
P_r=zeros(1,n);
for i=1:n   
%     if Rho(i)>=Rho_0
%         P_r(i)=0;
%     else
%         P_r(i)=A/(1+exp(-B*Rho(i)^2));
%     end
    P_r(i)=A/(1+exp(-B*Rho(i)^2));
end
end

