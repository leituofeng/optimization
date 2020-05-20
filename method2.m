clear;clc;close all;
rng('default');
m=500;n=1000;
a=randn(n,m);
b=sign(rand(m,1)-0.5);
w(1,:)=zeros(1,n);%假设初始的权值
c(1)=1;%初始的C
%第一次沿着负梯度方向搜索

%% %对于w的梯度
temp=0;
for i=1:m       
     temp=temp+(-b(i)*a(:,i)*exp(-b(i)*(w(1,:)*a(:,i)+c(1))))/(1+exp(-b(i)*(w(1,:)*a(:,i)+c(1))));
end
temp=temp./m;%初始点的梯度方向
g(:,1)=temp;
d(:,1)=temp;%初始点的搜索方向
buchang=0.2;
w(2,:)=w(1,:)-buchang*temp';
%% 对于C的梯度
temp0=0;
for i=1:m
    temp0=temp0+(-b(i)*exp(-b(i)*(w(1,:)*a(:,i)+c(1))))/(1+exp(-b(i)*(w(1,:)*a(:,i)+c(1))));
end
temp0=temp0/m;
c_d(1)=temp0;%c点处的梯度
c_g(1)=temp0;
c(2)=c(1)-buchang*c(1);
%% 开始迭代
k=2;
while(norm(temp)>=1e-2)
    temp=0;
    for i=1:m       
     temp=temp+(-b(i)*a(:,i)*exp(-b(i)*(w(k,:)*a(:,i)+c(1))))/(1+exp(-b(i)*(w(k,:)*a(:,i)+c(k))));
    end
    temp=temp/m;
    g(:,k)=temp;
    y(k-1)=norm(temp);
    belta=norm(g(:,k))/norm(g(:,k-1));
    d(:,k)=-temp+belta*d(:,k-1);
    w(k+1,:)=w(k,:)+buchang*d(:,k)';
    temp0=0;
    for i=1:m
        temp0=temp0+(-b(i)*exp(-b(i)*(w(k,:)*a(:,i)+c(k))))/(1+exp(-b(i)*(w(k,:)*a(:,i)+c(k))));
    end
    temp0=temp0/m;
    c_g(k)=temp0;
    c_belta=norm(temp0)/norm(c_g(k-1));
    c_d(k)=-temp0+c_belta*(c_d(k-1));
    c(k+1)=c(k)+buchang*c_d(k);
    k=k+1;
end
plot(y),xlabel('n'),ylabel('Euclidean norm');title('采用共轭梯度法')