function [scr_code mazhu]=scramble_seek(X,zhentou,raomazhu)
for i=1:8
    aaa=scramble(16*(8*raomazhu+i-1),1);
    m(i)=abs(X(zhentou:8:zhentou+2560*3*8-1)*aaa(1:2560*3)');  %��һ��֡���ݲ��ֽ�����أ�����š���ΪƵƫ�Ĵ��ڣ������˴󳤶ȵ�����ۼ�
end
ii=[0 1 2 3 4 5 6 7];
figure
plot(ii,m);xlim([-1 8]);title('������������8�������ƥ��ͼ');
scr_code=16*(8*raomazhu+find(m==max(m))-1);
mazhu=raomazhu;
%%%%%%%��ΪCPICH ʹ�õ���0��256�׵���Ƶ�룬����Ƶ��Ϊȫ1���У���������̼�Ϊÿ256������ӣ����䱾��ҲΪȫ1����.
%%%%%%%���ʹ��abs(X(zhentou:zhentou+38399)*scramble(raomazhu+(16*(i-1)),1))'һ��֡����ֱ�����Ϊͬ������ͨ���˷���ƥ�����ֵ��Ϊ��ȷ������