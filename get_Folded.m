function [F]=get_Folded(M,N)
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
F=M;
display('Folded');
printMatrix(F);
end