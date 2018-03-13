function data_jk=code_desp(data,ovsf_code)   %定义解扩函数，输入数据data， 时间帧头位置集合zhentou,时隙位置集合slot_ser，扩频码ovsf_code
N1=length(data);                          %解扰后的帧数
N2=length(ovsf_code);                %扩频码长度，阶数
data_tr=reshape(data,N2,N1/N2);    %按扩频阶数转置
data_jk=ovsf_code*data_tr;        %矩阵的乘法
% data_jk=zeros(1,2*length(data_jk1));
% data_jk(1:2:end)=real(data_jk1);
% data_jk(2:2:end)=imag(data_jk1);        %实部虚部组合奇偶e 
end