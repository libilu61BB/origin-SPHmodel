function P_r = P_repul(Rho,Rho_0)
%P_repul �����ų�������ÿ�����˵�ѹǿ
%   Rho �����ӵ��ܶ�
%   Rho_0 ���ӵ��ٽ��ܶ�
%   A ѹǿ����ϵ��
%   B ѹǿ����ϵ��
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

