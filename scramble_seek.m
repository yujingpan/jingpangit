function [scr_code mazhu]=scramble_seek(X,zhentou,raomazhu)
for i=1:8
    aaa=scramble(16*(8*raomazhu+i-1),1);
    m(i)=abs(X(zhentou:8:zhentou+2560*3*8-1)*aaa(1:2560*3)');  %对一个帧数据部分进行相关，扰码号。因为频偏的存在，不适宜大长度的相干累加
end
ii=[0 1 2 3 4 5 6 7];
figure
plot(ii,m);xlim([-1 8]);title('该主扰码组内8个扰码的匹配图');
scr_code=16*(8*raomazhu+find(m==max(m))-1);
mazhu=raomazhu;
%%%%%%%以为CPICH 使用的是0号256阶的扩频码，该扩频码为全1序列，其解扩过程即为每256序列相加，而其本身也为全1序列.
%%%%%%%因此使用abs(X(zhentou:zhentou+38399)*scramble(raomazhu+(16*(i-1)),1))'一个帧序列直接相加为同样处理，通过此方法匹配最大值即为正确的扰码