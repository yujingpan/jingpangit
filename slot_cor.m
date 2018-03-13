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

for i=1:N                                %N��ʱ϶      
    for a=1:16 
     threetoone=[abs(X(slot_ser(i)-1:slot_ser(i)+254)*ssc(a,:)'),abs(X(slot_ser(i):slot_ser(i)+255)*ssc(a,:)'),abs(X(slot_ser(i)+1:slot_ser(i)+256)*ssc(a,:)')];      %slot_ser(i)��Ϊ��i��ʱ϶����ʼλ��
     ppz(i,a)=max(threetoone);
    end
    ssc_code(i)=find(ppz(i,:)>=max(ppz(i,:)));                 %���ص�i���������λ�ã���a��ֵ��ΪSSC��ı��1��16֮��
    threetoone=[abs(X(slot_ser(i)-1:slot_ser(i)+254)*ssc(ssc_code(i),:)'),abs(X(slot_ser(i):slot_ser(i)+255)*ssc(ssc_code(i),:)'),abs(X(slot_ser(i)+1:slot_ser(i)+256)*ssc(ssc_code(i),:)')];
    pianchazhi(i)=(-2+find(threetoone==max(threetoone)));      %Ѱ��ʱ϶ͷ��ƫ��λ��
end
    pianchazhi;                                              %����ƫ��ֵ                                      
    slot_ser_xiuzheng=slot_ser+pianchazhi;                 %����������ʱ϶���