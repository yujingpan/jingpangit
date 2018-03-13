function c=scramble(n, len);
if(n>8191||n<0)
    error('WCDMA系统的扰码号为0-8191，请重新输入');
end

%返回扰码编号为n的扰码，长度为2帧
x=zeros(1,262142);
x(1,1)=1;
y=ones(1,262142);
z=zeros(1,262142);
c=zeros(1,38400*len);
N=38400*len;
%262142-18=262124;
L=262124;%2.^18-2
M=131072;%2.^17
 
for i=1:L;
    tempx=x(1,i+7)+x(1,i);
    x(1,i+18)=mod(tempx,2);
    tempy=y(1,i+10)+y(1,i+7)+y(1,i+5)+y(1,i);
    y(1,i+18)=mod(tempy,2);
end
for i=1:L
    if(i+n) >262124
        i1 = i + n - 262124;
    else
        i1 = i + n;
    end
    tempz=x(1,i1)+y(1,i); %第n号扰码
    z(1,i)=mod(tempz,2);
    if z(1,i)==1
        z(1,i)=-1;
    else 
        z(1,i)=1;
    end
end
for i=1:N;
    c(1,i)=z(1,i)+j*z(1,i+M);
end 
















