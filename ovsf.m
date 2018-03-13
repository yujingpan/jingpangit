function c=ovsf(sf_length) % sf_length为输入的扩频因子变量
c=zeros(sf_length,sf_length);
ctemp=zeros(sf_length,sf_length);
c(1,1)=1;
for i=1:log2(sf_length);
    itemp=2^i;
    ctemp=c;
    for j=1:itemp/2;
        c(2*j-1,1:itemp/2)=ctemp(j,1:itemp/2);
        c(2*j-1,itemp/2+1:itemp)=ctemp(j,1:itemp/2);
        c(2*j,1:itemp/2)=ctemp(j,1:itemp/2);
        c(2*j,itemp/2+1:itemp)=-ctemp(j,1:itemp/2);
    end
end
