function data_re=code_desc128_1(X,zhentou,scrcode)    
M=zhentou-1;
N=round(length(X)/(8*38400)-0.5)-1;
OV=128;
xixi=0;
xixi_xiu=zeros(1,38400*N/OV+1);
% xixi_xiu=0;
scrcode=kron(ones(1,N),scrcode);
data_re_1=zeros(1,38400*N);  
for i=1:38400*N/OV
cpe=[abs(X(OV*8*(i-1)+0+M+xixi_xiu(i):8:OV*8*i-8+M+xixi_xiu(i))*scrcode(OV*(i-1)+1:OV*i)'),abs(X(OV*8*(i-1)+1+M+xixi_xiu(i):8:OV*8*i-7+M+xixi_xiu(i))*scrcode(OV*(i-1)+1:OV*i)'),abs(X(OV*8*(i-1)+2+M+xixi_xiu(i):8:OV*i*8-6+M+xixi_xiu(i))*scrcode(OV*(i-1)+1:OV*i)')];
compe(i)=max(cpe);
xixi=1*(find(cpe==max(cpe))-2);
xixi_xiu(i+1)=xixi_xiu(i)+xixi;
data_re_1(OV*(i-1)+1:OV*i)=X(OV*8*(i-1)+1+M+xixi_xiu(i+1):8:OV*8*i+M-7+xixi_xiu(i+1)).*conj(scrcode(OV*(i-1)+1:OV*i));

end
figure;plot(xixi_xiu);
data_re= data_re_1;