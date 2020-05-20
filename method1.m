clear;clc;close all;
rng('default');
m=500;n=1000;
a=randn(n,m);
b=sign(rand(m,1)-0.5);
w(1,:)=zeros(1,n);%假设初始的权值
j=1;
c(j)=1;%假设初始的c
alpha=0.2;%步长
temp(:,j)=ones(1000,1);%初始梯度设大点
tidu(:,1)=ones(1000,1);
while(norm(tidu(:,j))>=1e-4)
    f(j)=0;
    for i=1:m
        f(j)=f(j)+log(1+exp(-b(i)*(w(j,:)*a(:,i)+c(j))));
    end
    f(j)=f(j)/m;
    temp(:,j)=0;
    for i=1:m       %对于w的梯度
        temp(:,j+1)=temp(:,j)+(-b(i)*a(:,i)*exp(-b(i)*(w(j,:)*a(:,i)+c(j))))/(1+exp(-b(i)*(w(j,:)*a(:,i)+c(j))));
    end
    tidu(:,j+1)=temp(:,j+1)./m;
    y(j)=norm(tidu(:,j+1));
    w(j+1,:)=w(j,:)-alpha*tidu(:,j+1)';
    temp0=0;
    
    for i=1:m
        temp0=temp0+(-b(i)*exp(-b(i)*(w(j,:)*a(:,i)+c(j))))/(1+exp(-b(i)*(w(j,:)*a(:,i)+c(j))));
    end
    temp0=temp0/m;
    c(j+1)=c(j)-alpha*c(j);
    j=j+1;
end
plot(y),xlabel('n'),ylabel('Euclidean norm');title('采用最速下降法')