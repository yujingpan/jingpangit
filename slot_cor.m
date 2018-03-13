function [slot_ser_xiuzheng,pianchazhi]=slot_cor(X,slot_ser)
N=length(slot_ser);
b=[1 1 1 1 1 1 -1 -1 -1 1 -1 1 -1 1 1 -1];
z=[b b b -b b b -b -b b -b b -b -b -b -b -b];
H=(1);
for k=1:8
    H=[H H;H -H];
end
H1=H(1:16:241,:);
ssc=zeros(16,256);

for k=1:16
    ssc(k,:)=(1+sqrt(-1))*H1(k,:).*z;
end

for i=1:N                                %N个时隙      
    for a=1:16 
     threetoone=[abs(X(slot_ser(i)-1:slot_ser(i)+254)*ssc(a,:)'),abs(X(slot_ser(i):slot_ser(i)+255)*ssc(a,:)'),abs(X(slot_ser(i)+1:slot_ser(i)+256)*ssc(a,:)')];      %slot_ser(i)作为第i个时隙的起始位置
     ppz(i,a)=max(threetoone);
    end
    ssc_code(i)=find(ppz(i,:)>=max(ppz(i,:)));                 %返回第i行最大数的位置，即a点值，为SSC码的编号1：16之间
    threetoone=[abs(X(slot_ser(i)-1:slot_ser(i)+254)*ssc(ssc_code(i),:)'),abs(X(slot_ser(i):slot_ser(i)+255)*ssc(ssc_code(i),:)'),abs(X(slot_ser(i)+1:slot_ser(i)+256)*ssc(ssc_code(i),:)')];
    pianchazhi(i)=(-2+find(threetoone==max(threetoone)));      %寻找时隙头的偏差位置
end
    pianchazhi;                                              %返回偏差值                                      
    slot_ser_xiuzheng=slot_ser+pianchazhi;                 %返回修正后时隙起点