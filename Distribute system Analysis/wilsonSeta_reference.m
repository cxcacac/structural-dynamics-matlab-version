function [ua,ua1,ua2]=wilson(m,k,c,pt,u,u1,u2,n,theta,dertat) 
a0=6/((theta*dertat)^2); 
a1=3/(theta*dertat); 
a2=2*a1; 
a3=theta*dertat/2; 
a4=a0/theta; 
a5=-a2/theta; 
a6=1-3/theta; 
a7=dertat/2; 
a8=dertat^2/6; 
k=k+a0*m+a1*c; 
[r,s]=ldl(k); 

for i=0:n
    pta=pt+m*(a0*u+a2*u1+2*u2)+c*(a1*u+2*u1+a3*u2); 
    uaa=r\pta; 
    uaa=s\uaa; 
    uaa=r'\uaa; 
    ua2(:,i+1)=a4*(uaa-u)+a5*u1+a6*u2; 
    ua1(:,i+1)=u1+a7*(ua2(:,i+1)+u2); 
    ua(:,i+1)=u+dertat*u1+a8*(ua2(:,i+1)+2*u2); 
    u=ua(:,i+1); 
    u1=ua1(:,i+1); 
    u2=ua2(:,i+1); 
end