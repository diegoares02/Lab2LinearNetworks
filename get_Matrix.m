function [M]=get_Matrix(A,B,C,D,Ps,epsilon,N)
Y11_num=D;
Y12_num= (-1i*Ps/epsilon);
Y21_num= (-1i*Ps/epsilon);
Y22_num=A;
Y_den=B;

[r11,lambda,K]=residue(Y11_num,Y_den);
[r12,lambda,K]=residue(Y12_num,Y_den);
[r21,lambda,K]=residue(Y21_num,Y_den);
[r22,lambda,K]=residue(Y22_num,Y_den);

if isempty(K)
    K_inf=0;
end
lambda=real(-1i*lambda);
Msl=K_inf;
Ck=1;
Bk=-lambda;
T1k=1i*r21./sqrt(r22);
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
display('Initial Matrix');
printMatrix(M);
end