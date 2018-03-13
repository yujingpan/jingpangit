%Lowpass filter 设计一个低通滤波器
function y=lowp1(x,bw)
fs =15000;
fchip = 3840000;    % 3.84 Mhz
%bw =  4000000;      % 带宽   Mhz
% fc = 20000000;  
    Ap= 0.5;
    Ast= 20;
    BW = bw;
    Fc1 =bw;
    Fs = fs;
    
	fp= Fc1 - 0.05*BW/2;
    fst= Fc1 + 0.05*BW/2;
    ap= Ap;
    ast= Ast;
    hs1 = fdesign.lowpass( 'Fp,Fst,Ap,Ast',fp,fst,ap,ast,Fs);
    bpType = 'cheby1';
    bpf1 = design(hs1, bpType);
    y=filter(bpf1,x);  
end