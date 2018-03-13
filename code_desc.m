function data_re=code_desc(X,zhentou,scrcode)    
M=zhentou-8;
N=round(length(X)/(38400*8)-0.5)-1;
OV=256;
xixi=0;
xixi_xiu=0;
scrcode=kron(ones(1,N),scrcode);
data_re_1=zeros(1,38400*N);  
for i=1:38400*N/OV
cpe=[abs(X(OV*8*(i-1)+0+M+xixi_xiu:8:OV*8*i-8+M+xixi_xiu)*scrcode(OV*(i-1)+1:OV*i)'),abs(X(OV*8*(i-1)+8+M+xixi_xiu:8:OV*8*i+M+xixi_xiu)*scrcode(OV*(i-1)+1:OV*i)'),abs(X(OV*8*(i-1)+16+M+xixi_xiu:8:OV*i*8+8+M+xixi_xiu)*scrcode(OV*(i-1)+1:OV*i)')];
compe(i)=max(cpe);
xixi=8*(find(cpe==max(cpe))-2);
xixi_xiu=xixi_xiu+xixi;
data_re_1(OV*(i-1)+1:OV*i)=X(OV*8*(i-1)+8+M+xixi_xiu:8:OV*8*i+M+xixi_xiu).*conj(scrcode(OV*(i-1)+1:OV*i));
end
data_re= data_re_1;