function data_jk=code_desp(data,ovsf_code)   %���������������������data�� ʱ��֡ͷλ�ü���zhentou,ʱ϶λ�ü���slot_ser����Ƶ��ovsf_code
N1=length(data);                          %���ź��֡��
N2=length(ovsf_code);                %��Ƶ�볤�ȣ�����
data_tr=reshape(data,N2,N1/N2);    %����Ƶ����ת��
data_jk=ovsf_code*data_tr;        %����ĳ˷�
% data_jk=zeros(1,2*length(data_jk1));
% data_jk(1:2:end)=real(data_jk1);
% data_jk(2:2:end)=imag(data_jk1);        %ʵ���鲿�����że 
end