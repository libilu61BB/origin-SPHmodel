function P_e = P_extru(Rho,Rho_0)
%P_extru ���������໥�Ӵ���ѹ������ѹǿ
%   Rho �����ӵ��ܶ�
%   Rho_0 ���ӵ��ٽ��ܶ�
%   Radius ���������ӵİ뾶
n=length(Rho);
K=1000; %����
P_e=zeros(1,n);
for i=1:n
    P_e(i)=K*(Rho(i)-Rho_0)/Rho_0;
    if P_e(i)<0
        P_e(i)=0;%��������Ϊ������ѹǿ����Ϊ0
    end
end
end