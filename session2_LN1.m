clear all
epsilon=43.4176;
RL=20;
P_s=[1i 0 4.85*1i 0 5.6644*1i];
F_s=[1 0 2.1412 0 1.4691 0 0.3369 0 0.0132];
E_s=[1 1.9462 4.035 4.718 4.8027 3.4133 1.8758 0.6764 0.1322];
N=length(E_s)-1; %order of the filter 

%ABCD matrix
for k=1:N+1
    if  mod(k,2)~=0
        A_s(N+2-k)=1i*imag(E_s(N+2-k)+F_s(N+2-k));
        B_s(N+2-k)=real(E_s(N+2-k)+F_s(N+2-k));
        C_s(N+2-k)=real(E_s(N+2-k)-F_s(N+2-k));
        D_s(N+2-k)=1i*imag(E_s(N+2-k)-F_s(N+2-k));
    end
    if  mod(k,2)==0   
        A_s(N+2-k)=real(E_s(N+2-k)+F_s(N+2-k));
        B_s(N+2-k)=1i*imag(E_s(N+2-k)+F_s(N+2-k));
        C_s(N+2-k)=1i*imag(E_s(N+2-k)-F_s(N+2-k));
        D_s(N+2-k)=real(E_s(N+2-k)-F_s(N+2-k));        
    end
    
end

%Admittance matrix
Y11_num=D_s;
Y12_num= -(1i*P_s)/epsilon;
Y21_num= -(1i*P_s)/epsilon;
Y22_num=A_s;
Y_den=B_s;

[r11,lambda,K]=residue(Y11_num,Y_den);
[r12,lambda,K]=residue(Y12_num,Y_den);
[r21,lambda,K]=residue(Y21_num,Y_den);
[r22,lambda,K]=residue(Y22_num,Y_den);

if isempty(K)
    K_inf=0;
end

Msl=K_inf;
Ck=1;
Bk=-lambda;
T1k=r21./sqrt(r22);
Tnk=sqrt(r22);
M=zeros(N+2);

for k=2:N+1
         M(k,k)=Bk(k-1);
         M(k,1)=T1k(k-1);
         M(1,k)=T1k(k-1);
         M(N+2,k)=Tnk(k-1);
         M(k,N+2)=Tnk(k-1);   
end

M(1,N+2)=Msl;
M(N+2,1)=Msl;

%elements to eliminate= s2,s3,s4,s5,L2,L4,L3,13,14,53
R=zeros(N+2);
%for k=1:29
 %   M(:,k)=R(:,k)*M(:,k-1)*R(:,k).';













