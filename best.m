function ww=best(X,zhentou,scrcode)    
M=zhentou-4;
for i=1:8
   m(i)=abs(X(M+i:8:M+i+2560*4*8-1)*scrcode(1:2560*4)'); 
end
% plot(m);title('最佳取样位置')
ww=find(m==max(m));
