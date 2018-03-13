clear all
fs=18.75e6;                               %������
fc=fs/4-0.4e6-1130;                       %�ز�Ƶ��
pname = 'E:\�����ļ�\';
fname = 'WCDMA 2143Mhz_S18.75MHZ_B8M_16bit.dat';
fid=fopen([pname fname],'r');             %    �ļ� ָ��
nsamp=8;                                  %�ز�������
nFrame=5;                                 % ֡��
nslot=15*nFrame;                          %ʱ϶��
Lslot1=2560*5;                            %����fs/3.84e6=5 ��ʱ϶����
L_slot=2560*nsamp;                        % �ز�����ʱ϶����
N=Lslot1*nslot;                           % ȡ�� 
N_re=1572864;                             %ȡ��N �ϲ����� 8��ʱ��ó���
N_f=12;                                   %Сѭ������
s=zeros(1,N_f*N_re);        
N_b=120;                                    % ��ѭ������   ��֡��=nFrame*N_f*N_b        
%%
%������ȡ��N_f�Σ��ϲ�
for jj=1:N_b                                   % �ļ�̫����Ҫ�ִ���������
    N_SE(1)=0;                                 %ÿ��ȡ��ƫ��λ��
    nj(1)=0;                                   %������ǰ��ȡ��֡��  
for i=1:N_f                                         %N_f��ѭ����һ�ζ�ȡnFrame*N_f֡����
fseek(fid,2*N*(i-1)+N_SE(jj),'bof');                 % �ƶ���ʼλ��
data_be=fread(fid,N,'int16');                     % ָ���ļ�β
n=0:N-1;t=(n/fs); 
f1=n*fs/N;                    
dx=1*pi/8;                                    %��ƫ��
data_be=(hilbert(data_be).').*exp(-sqrt(-1)*(fc*2*pi*t+dx+(N_SE(jj)/fs/2+N*(i-1)/fs)*fc*2*pi));             % ϣ�����ر任    �±�Ƶ
data_be=lowp(data_be,2600000);                %��ͨ�˲�                     
data_1=resample(data_be,3.84e6*nsamp,fs);     %�ز�����������3.84e6*nsamp
s(N_re*(i-1)+1:N_re*i)=data_1;      
end        
clear data_be;
clear data_1

%%
%%%%%%%%%%%%%%%%%%  PSC ͬ�� %%%%%%%%%%%%%%%%%%
%  u=[1,1,1,1,1,1,-1,-1,1,-1,1,-1,1,-1,-1,1];
%  v=[1,1,1,-1,-1,1,-1,-1,1,1,1,-1,1,-1,1,1];
% 
% PSC=((1+sqrt(-1))*(kron(v,u)));    
% PSC=kron(PSC,[1 zeros(1,nsamp-1) ]);
% for i=1:L_slot*30;
%     corr(i)=abs(s(i:i+256*(nsamp/2)-1)*(PSC(1:256*4)'))+abs(s(i+256*(nsamp/2):i+256*nsamp-1)*(PSC(1+256*4:256*8)'));
% end
% corr_1= reshape(abs(corr),2560*8,30);
% corr_2= sum(corr_1,2);
% figure;
% subplot(211)
% plot(corr);  title('ʱ϶���ͼ')
% subplot(212)
% plot(corr_2);  title('ʱ϶�ۼ�ͼ')
% slot_start=find(corr_2==max(corr_2));     %ȷ����һ��ʱ϶λ��3414��ʵ���ϻ�������һ���Ƚϴ�ĵ�2107�� ����һ�����ʽϴ��С����ʱ϶ͷλ��  %                                            
slot_start=3414;
                   %%%slot_start���ѡ�񣬾�����ѡ��ȡ�ĸ�С������Ϣ

%%
%%%%%%%%%%%�ڶ���          
% b=[1 1 1 1 1 1 -1 -1 -1 1 -1 1 -1 1 1 -1];
% Z=[b b b -b b b -b -b b -b b -b -b -b -b -b];
% H0 = 1;
% H1 = hardmard_gen(H0);
% H2 = hardmard_gen(H1);
% H3 = hardmard_gen(H2);
% H4 = hardmard_gen(H3);
% H5 = hardmard_gen(H4);
% H6 = hardmard_gen(H5);
% H7 = hardmard_gen(H6);
% H8 = hardmard_gen(H7); 
% hardmardSequence = H8(1:16:256,:);
% Cssc = kron(ones(16,256),(1+1i)).*(kron(ones(16,1),Z).*hardmardSequence);
% % 
% % % FHT  ������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% ZZ = (1-1i)*Z;
% rFHT = zeros(15,16);
%  
% for j = 0:2          %%%%%%%�����ۼӵ�֡����
%     frame_start = slot_start + j*L_slot*15;
%     for i = 0:14
%         slot_start_now = frame_start + i*L_slot;
%         
%         ss = s(slot_start_now:nsamp:slot_start_now+nsamp*length(ZZ)-1).*ZZ;   
%         
%         sss = reshape(ss,16,length(ss)/16);
%         sss = sum(sss,1);
%         rFHT(i+1,:) = rFHT(i+1,:) + abs((H4*sss.')/length(sss)).';          %16��SSCƥ��
%     end
% end
% rFHT = abs(rFHT);
% [maxb,max_p] = max(rFHT,[],2);         %�ɷ��ؾ���������������ο�
% Cssc_pos = max_p.';
% Cssc_pos_matrix = zeros(length(Cssc_pos));
% for i = 1:15
%     for j = 1:15
%         Cssc_pos_matrix(i,j) = Cssc_pos(mod(i+j-2,15)+1);
%     end
% end
%    
% sscTable = [
% 1	1	2	8	9	10	15	8	10	16	2	7	15	7	16
% 1	1	5	16	7	3	14	16	3	10	5	12	14	12	10
% 1	2	1	15	5	5	12	16	6	11	2	16	11	15	12
% 1	2	3	1	8	6	5	2	5	8	4	4	6	3	7
% 1	2	16	6	6	11	15	5	12	1	15	12	16	11	2
% 1	3	4	7	4	1	5	5	3	6	2	8	7	6	8
% 1	4	11	3	4	10	9	2	11	2	10	12	12	9	3
% 1	5	6	6	14	9	10	2	13	9	2	5	14	1	13
% 1	6	10	10	4	11	7	13	16	11	13	6	4	1	16
% 1	6	13	2	14	2	6	5	5	13	10	9	1	14	10
% 1	7	8	5	7	2	4	3	8	3	2	6	6	4	5
% 1	7	10	9	16	7	9	15	1	8	16	8	15	2	2
% 1	8	12	9	9	4	13	16	5	1	13	5	12	4	8
% 1	8	14	10	14	1	15	15	8	5	11	4	10	5	4
% 1	9	2	15	15	16	10	7	8	1	10	8	2	16	9
% 1	9	15	6	16	2	13	14	10	11	7	4	5	12	3
% 1	10	9	11	15	7	6	4	16	5	2	12	13	3	14
% 1	11	14	4	13	2	9	10	12	16	8	5	3	15	6
% 1	12	12	13	14	7	2	8	14	2	1	13	11	8	11
% 1	12	15	5	4	14	3	16	7	8	6	2	10	11	13
% 1	15	4	3	7	6	10	13	12	5	14	16	8	2	11
% 1	16	3	12	11	9	13	5	8	2	14	7	4	10	15
% 2	2	5	10	16	11	3	10	11	8	5	13	3	13	8
% 2	2	12	3	15	5	8	3	5	14	12	9	8	9	14
% 2	3	6	16	12	16	3	13	13	6	7	9	2	12	7
% 2	3	8	2	9	15	14	3	14	9	5	5	15	8	12
% 2	4	7	9	5	4	9	11	2	14	5	14	11	16	16
% 2	4	13	12	12	7	15	10	5	2	15	5	13	7	4
% 2	5	9	9	3	12	8	14	15	12	14	5	3	2	15
% 2	5	11	7	2	11	9	4	16	7	16	9	14	14	4
% 2	6	2	13	3	3	12	9	7	16	6	9	16	13	12
% 2	6	9	7	7	16	13	3	12	2	13	12	9	16	6
% 2	7	12	15	2	12	4	10	13	15	13	4	5	5	10
% 2	7	14	16	5	9	2	9	16	11	11	5	7	4	14
% 2	8	5	12	5	2	14	14	8	15	3	9	12	15	9
% 2	9	13	4	2	13	8	11	6	4	6	8	15	15	11
% 2	10	3	2	13	16	8	10	8	13	11	11	16	3	5
% 2	11	15	3	11	6	14	10	15	10	6	7	7	14	3
% 2	16	4	5	16	14	7	11	4	11	14	9	9	7	5
% 3	3	4	6	11	12	13	6	12	14	4	5	13	5	14
% 3	3	6	5	16	9	15	5	9	10	6	4	15	4	10
% 3	4	5	14	4	6	12	13	5	13	6	11	11	12	14
% 3	4	9	16	10	4	16	15	3	5	10	5	15	6	6
% 3	4	16	10	5	10	4	9	9	16	15	6	3	5	15
% 3	5	12	11	14	5	11	13	3	6	14	6	13	4	4
% 3	6	4	10	6	5	9	15	4	15	5	16	16	9	10
% 3	7	8	8	16	11	12	4	15	11	4	7	16	3	15
% 3	7	16	11	4	15	3	15	11	12	12	4	7	8	16
% 3	8	7	15	4	8	15	12	3	16	4	16	12	11	11
% 3	8	15	4	16	4	8	7	7	15	12	11	3	16	12
% 3	10	10	15	16	5	4	6	16	4	3	15	9	6	9
% 3	13	11	5	4	12	4	11	6	6	5	3	14	13	12
% 3	14	7	9	14	10	13	8	7	8	10	4	4	13	9
% 5	5	8	14	16	13	6	14	13	7	8	15	6	15	7
% 5	6	11	7	10	8	5	8	7	12	12	10	6	9	11
% 5	6	13	8	13	5	7	7	6	16	14	15	8	16	15
% 5	7	9	10	7	11	6	12	9	12	11	8	8	6	10
% 5	9	6	8	10	9	8	12	5	11	10	11	12	7	7
% 5	10	10	12	8	11	9	7	8	9	5	12	6	7	6
% 5	10	12	6	5	12	8	9	7	6	7	8	11	11	9
% 5	13	15	15	14	8	6	7	16	8	7	13	14	5	16
% 9	10	13	10	11	15	15	9	16	12	14	13	16	14	11
% 9	11	12	15	12	9	13	13	11	14	10	16	15	14	16
% 9	12	10	15	13	14	9	14	15	11	11	13	12	16	10];
% 
% sscTable_match_pos = 0;                                    % ��ʼ�����������
% sscTable_match_num = 0;                                    
% frame_start = 0;                                           % ��ʼ��֡��ʼλ��                  
% for i = 1:15
%     for j = 1:64
%         match_num = 0;
%         for k=1:15
%             if(Cssc_pos_matrix(i,k) == sscTable(j,k))
%                 match_num = match_num+1;
%             end
%         end
%         if(sscTable_match_num<match_num)
%             sscTable_match_pos = j;
%             frame_start = (i-1)*L_slot;
%             sscTable_match_num = match_num;
%         end
%     end
% end

% frame_location(1)=frame_start+slot_start;          %֡��ʼλ��
% sscTable_match_pos = sscTable_match_pos - 1;   % ����������0��ʼ
frame_location(1)=269654;          
% sscTable_match_pos_use(1)= sscTable_match_pos;      %%% ����һ��С����slot_start=2107�� frame_start =204506;sscTable_match_pos = 10;

%%
%%%%%%%%%������ ������ƥ��
% [scr_code sscTable_match_pos_use(jj+1)] =scramble_seek(s,frame_location(jj),sscTable_match_pos_use(jj));   %%%%�������4336  �� ��һ��С����1296��   
% ww=best(s,frame_location(jj),scramble(scr_code,1));
scr_code=4336;           
%%
%%%����
[data_des,frame_location(jj+1),nj(jj)]=code_desc1(s,frame_location(jj),scramble(scr_code,1));
%��һ��ȡ��֡λ��         frame_location(jj+1)
% ��һ��ȡ����ǰ��ȡ֡��  nj(jj)

%%
%�Ե�Ƶ�ŵ�P-CPICH����
data_cpi=reshape(data_des,256,length(data_des)/256);
aaa=sum(data_cpi,1);        % ��Ƶ��Ƶ�� C256,0
daopin1(9000*(jj-1)+1:9000*jj)=aaa;
esti_co=(1+1i)./aaa;      % �ŵ����Ʋ�������
% figure;
% subplot(311);plot(real(aaa));title('��Ƶʵ��');
% subplot(312);plot(imag(aaa));title('��Ƶ�鲿');
% fs2=15000;n=0:length(aaa)-1;t1=(n/fs2); f2=n*fs2/length(aaa);   
% subplot(313);plot(f2-fs2/2,fftshift(abs(fft(real(aaa)))));xlabel('HZ');title('Ƶƫ')
% figure;
% plot(abs(aaa));
%%
%�������ŵ����� ����BCH��P-CCPCH�ŵ�
ssss=ovsf(256);
OVSF=ssss(2,:);     % BCH ����Ƶ�� C256,1
s_channel_1 = reshape(data_des,256,length(data_des)/256).';
s_channel_1 = s_channel_1*OVSF';
s_channel_1 = s_channel_1.*(esti_co).';    
data_xx=reshape(s_channel_1,10,length(s_channel_1)/10);
data_yy=data_xx(2:10,:);
data_yy=reshape(data_yy,1,length(s_channel_1)*9/10);
BCH=ones(1,length(s_channel_1)*18/10);
BCH(1:2:end)=real(data_yy);
BCH(2:2:end)=imag(data_yy);
BCH=sign(BCH);
BCH=(-BCH+1)/2;
data_w1(length(BCH)*(jj-1)+1:length(BCH)*jj)=BCH;

%%
%�������ŵ�����   ����Ѱ��ָʾPICH��S-CCPCH�ŵ� 
ssss=ovsf(256);
OVSF_2=ssss(4,:);       %   PICH����Ƶ�� C256,3
s_channel_2 = reshape(data_des,256,length(data_des)/256).';
s_channel_2 = s_channel_2*OVSF_2';
s_channel_2 = s_channel_2.*(esti_co).';  
data_mi=ones(1,2*length(s_channel_2));data_mi(1:2:end)=real(s_channel_2);data_mi(2:2:end)=imag(s_channel_2);
data_pi(18000*(jj-1)+1:18000*jj)=data_mi;

%%
%�������ŵ�����  ���� PCH ��S-CCPCH�ŵ�      %��ʱû�з��������
ssss_2=ovsf(128);
OVSF_3 = ssss_2(5,:);            % PCH ��Ƶ�� C128,4
s_channel_3 = reshape(data_des,128,length(data_des)/128).';
s_channel_3 = s_channel_3*OVSF_3';
esti_co=resample(esti_co,2,1);
s_channel_3 = s_channel_3.*(esti_co).';
xunhu(18000*(jj-1)+1:18000*jj)=s_channel_3.';

clear esti_co��
N_SE(jj+1)=N_SE(jj)+2*N*N_f-187500*2-nj(jj)*2*187500;                       %��һ��ȡ����λָ����Ҫ�ƶ���λ��

end
save('data_cpich1','daopin1');
save('data_bch2','data_w1');          % ��������ŵ�
save('data_pich','data_pi'); 
save('data_pch','xunhu');          % ��������ŵ�





figure;subplot(211);plot([1:length(daopin1)]/150,real(daopin1)); title('real data of P-CPICH')
     subplot(212);plot([1:length(daopin1)]/150,imag(daopin1));xlabel('֡��');title('imag data of P-CPICH')
figure;plot([1:length(data_pi)]/300,data_pi); xlabel('֡��');title('data of S-CCPCH(PICH)')

figure;subplot(211);plot([1:length(xunhu)]/300,real(xunhu)); title('real data of P-CPICH')
     subplot(212);plot([1:length(xunhu)]/300,imag(xunhu));xlabel('֡��');title('imag data of P-CPICH')
