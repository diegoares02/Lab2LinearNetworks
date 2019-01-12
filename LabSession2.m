clear;clc;

epsilon=43.4176;
RL=20;
Ps=[1i 0 4.85*1i 0 5.6644*1i];
Fs=[1 0 2.1412 0 1.4691 0 0.3369 0 0.0132];
Es=[1 1.9462 4.035 4.718 4.8027 3.4133 1.8758 0.6764 0.1322];
N=length(Es)-1; %order of the filter 
%Admittance matrix
[A,B,C,D]=get_ABCD(Es,Fs);
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
original=M;
disp('Matrix M')
for k=1:N+2
    for j=1:N+2
        fprintf('%10.4f ',M(k,j));
    end
    fprintf('\n')
end

M(1,N+2)=Msl;
M(N+2,1)=Msl;
count=1;
display('    pivot  k  l  m  n  c   theta');
k=1;
c=-1;
 for j=N+1:-1:k+2
     R=eye(N+2);
     theta_r(j)=c*atan(M(k,j)/M(k,j-1));
     R(j-1,j-1)=cos(theta_r(j));
     R(j,j)=cos(theta_r(j));
     R(j,j-1)=sin(theta_r(j));
     R(j-1,j)=-sin(theta_r(j));
     M=R*M*R';
     fprintf('%d  [%d,%d]  %d  %d  %d  %d  %d   %f\n',count,j-1,j,k,j-1,k,j,c,theta_r(j))
     count=count+1;
 end
 
c=1;
k=10;
for j=3:8
    R=eye(N+2);
    theta_r(j)=c*atan(M(j,k)/M(j+1,k));
    R(j,j)=cos(theta_r(j));
    R(j+1,j+1)=cos(theta_r(j));
    R(j,j+1)=-sin(theta_r(j));
    R(j+1,j)=sin(theta_r(j));
    M=R*M*R';
    fprintf('%d  [%d,%d]  %d  %d  %d  %d  %d   %f\n',count,j,j+1,j,k,j+1,k,c,theta_r(j))
    count=count+1;
end


 k=2;
 c=-1;
for j=N:-1:k+2
    R=eye(N+2);
    theta_r(j)=c*atan(M(k,j)/M(k,j-1));
    R(j-1,j-1)=cos(theta_r(j));
    R(j,j)=cos(theta_r(j));
    R(j,j-1)=sin(theta_r(j));
    R(j-1,j)=-sin(theta_r(j));
    M=R*M*R';
    fprintf('%d  [%d,%d]  %d  %d  %d  %d  %d   %f\n',count,j-1,j,k,j-1,k,j,c,theta_r(j))
    count=count+1;
end
  
k=9; 
c=1;
for j=4:7
    R=eye(N+2);
    theta_r(j)=c*atan(M(j,k)/M(j+1,k));
    R(j+1,j+1)=cos(theta_r(j));
    R(j,j)=cos(theta_r(j));
    R(j,j+1)=-sin(theta_r(j));
    R(j+1,j)=sin(theta_r(j));
    M=R*M*R';
    fprintf('%d  [%d,%d]  %d  %d  %d  %d  %d   %f\n',count,j,j+1,j,k,j+1,k,c,theta_r(j))
    count=count+1;
end
 
k=3;
c=-1;
for j=N-1:-1:k+2
    R=eye(N+2);
    theta_r(j)=c*atan(M(k,j)/M(k,j-1));
    R(j-1,j-1)=cos(theta_r(j));
    R(j,j)=cos(theta_r(j));
    R(j,j-1)=sin(theta_r(j));
    R(j-1,j)=-sin(theta_r(j));
    M=R*M*R';
    fprintf('%d  [%d,%d]  %d  %d  %d  %d  %d   %f\n',count,j-1,j,k,j-1,k,j,c,theta_r(j))
    count=count+1;
end
   
k=8;
c=1;
for j=5:6
    R=eye(N+2);
    theta_r(j)=c*atan(M(j,k)/M(j+1,k));
    R(j+1,j+1)=cos(theta_r(j));
    R(j,j)=cos(theta_r(j));
    R(j,j+1)=-sin(theta_r(j));
    R(j+1,j)=sin(theta_r(j));
    M=R*M*R';
    fprintf('%d  [%d,%d]  %d  %d  %d  %d  %d   %f\n',count,j,j+1,j,k,j+1,k,c,theta_r(j))
    count=count+1;
end
k=4;
c=-1;
for j= N-2:-1:k+2
    theta_r=c*atan(M(k,j)/M(k,j-1));
    R=eye(N+2);
    R(6,6)=cos(theta_r);
    R(5,5)=cos(theta_r);
    R(6,5)=sin(theta_r);
    R(5,6)=-sin(theta_r);
    M=R*M*R';
    fprintf('%d  [%d,%d]  %d  %d  %d  %d  %d   %f\n',count,j-1,j,k,j-1,k,j,c,theta_r)
    count=count+1;
end
display('Folded');
for k=1:N+2
fprintf('%10.2f  ',M(k,:));
fprintf('\n');
end
folded=M;
a=M(3,8)*M(4,5)*M(5,6)-M(3,4)*M(6,7)*M(7,8)+M(3,8)*M(4,7)*M(6,7);
b=M(3,4)*M(4,7)*M(7,8)-M(3,8)*(M(4,5)^2-M(5,6)^2-M(6,7)^2+M(4,7)^2);
c=-M(3,8)*(M(4,7)*M(6,7)+M(4,5)*M(5,6));
t1=roots([a b c]);
Tk=[4 5 3];
Tl=[7 8 6];
Tm=[4 5 5];
Tn=[5 6 6];
Tc=[-1 -1 1];
i=[4 5 6 3];
j=[6 7 8 5];
theta_t=atan(t1(1));
CQ=folded;
R=eye(N+2);
R(i(1),i(1))=cos(theta_t);
R(j(1),j(1))=cos(theta_t);
R(j(1),i(1))=sin(theta_t);
R(i(1),j(1))=-sin(theta_t);
CQ=R*CQ*R';
for in=2:length(j)
    theta_t=Tc(in-1)*atan(M(Tk(in-1),Tl(in-1))/M(Tm(in-1),Tn(in-1)));
    R=eye(N+2);
    R(i(1),i(1))=cos(theta_t);
    R(j(1),j(1))=cos(theta_t);
    R(j(1),i(1))=sin(theta_t);
    R(i(1),j(1))=-sin(theta_t);
    CQ=R*CQ*R';
end
fprintf('Quartet Matrix\n');
for k=1:N+2
fprintf('%10.2f  ',CQ(k,:));
fprintf('\n');
end


W=eye(N+2);
W(1,1)=0;
W(N+2,N+2)=0;
R2=zeros(N+2);
R2(1,1)=1;
R2(N+2,N+2)=1;
indice=1;
for w=-5:0.01:5
    Am1=-1i*R2+w*W+M;
    Am2=-1i*R2+w*W+folded;
    Am3=-1i*R2+w*W+CQ;
    Am1=inv(Am1);
    Am2=inv(Am2);
    Am3=inv(Am3);
    S21t(indice)=2*1i*(Am1(N+2,1));
    S11t(indice)=-(1+2*1i*(Am1(1,1)));
    S21f(indice)=2*1i*(Am2(N+2,1));
    S11f(indice)=-(1+2*1i*(Am2(1,1)));
    S21q(indice)=2*1i*(Am3(N+2,1));
    S11q(indice)=-(1+2*1i*Am3(1,1));
    indice=indice+1;
end
d=-5:0.01:5;
S21tran=20*log10(abs(S21t));
S11tran=20*log10(abs(S11t));
S21fol=20*log10(abs(S21f));
S11fol=20*log10(abs(S11f));
S21qua=20*log10(abs(S21q));
S11qua=20*log10(abs(S11q));

figure(1)
plot(d,S21tran);
hold on
plot(d,S11tran,'r');
axis([min(d) max(d) -350 10])
title('Transversal Configuration')
xlabel('Frequency');
ylabel('Level(dB)');
grid on

figure(2)
plot(d,S21fol);
hold on
plot(d,S11fol,'r');
axis([min(d) max(d) -350 10])
title('Folded Configuration')
xlabel('Frequency');
ylabel('Level(dB)');
grid on

figure(3)
plot(d,S21qua);
hold on
plot(d,S11qua,'r');
axis([min(d) max(d) -350 10])
title('Quarter Configuration')
xlabel('Frequency');
ylabel('Level(dB)');
grid on