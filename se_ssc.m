function [code_group_number,frame_boder,phase_position]=se_ssc(x,slot_boder)
b=[1 1 1 1 1 1 -1 -1 -1 1 -1 1 -1 1 1 -1];
z=[b b b -b b b -b -b b -b b -b -b -b -b -b]; 
H=(1);
for k=1:8
    H=[H H;H -H];
end
H1=H(1:16:241,:);
ssc=zeros(16,256);
for k=1:16
    ssc(k,:)=(1+j)*H1(k,:).*z;
end
for k=1:15
    for n=1:16
        temp=slot_boder(k);
        correlation_value_16_15(n,k)= sum(ssc(n,:).*x(temp:(temp+255)));
    end
end


%次级同步码的分配表；
correlation_value_64_15=zeros(64,15);
table=[ 1	1	2	8	9	10	15	8	10	16	2	7	15	7	16; 
        1	1	5	16	7	3	14	16	3	10	5	12	14	12	10;
        1	2	1	15	5	5	12	16	6	11	2	16	11	15	12;
        1	2	3	1	8	6	5	2	5	8	4	4	6	3	7;
        1	2	16	6	6	11	15	5	12	1	15	12	16	11	2;
        1	3	4	7	4	1	5	5	3	6	2	8	7	6	8;
        1	4	11	3	4	10	9	2	11	2	10	12	12	9	3;
        1	5	6	6	14	9	10	2	13	9	2	5	14	1	13;
        1	6	10	10	4	11	7	13	16	11	13	6	4	1	16;
        1	6	13	2	14	2	6	5	5	13	10	9	1	14	10;
        1	7	8	5	7	2	4	3	8	3	2	6	6	4	5;
        1	7	10	9	16	7	9	15	1	8	16	8	15	2	2;
        1	8	12	9	9	4	13	16	5	1	13	5	12	4	8;
        1	8	14	10	14	1	15	15	8	5	11	4	10	5	4;
        1	9	2	15	15	16	10	7	8	1	10	8	2	16	9;
        1	9	15	6	16	2	13	14	10	11	7	4	5	12	3;
        1	10	9	11	15	7	6	4	16	5	2	12	13	3	14;
        1	11	14	4	13	2	9	10	12	16	8	5	3	15	6;
        1	12	12	13	14	7	2	8	14	2	1	13	11	8	11;
        1	12	15	5	4	14	3	16	7	8	6	2	10	11	13;
        1	15	4	3	7	6	10	13	12	5	14	16	8	2	11;
        1	16	3	12	11	9	13	5	8	2	14	7	4	10	15;
        2	2	5	10	16	11	3	10	11	8	5	13	3	13	8;
        2	2	12	3	15	5	8	3	5	14	12	9	8	9	14;
        2	3	6	16	12	16	3	13	13	6	7	9	2	12	7;
        2	3	8	2	9	15	14	3	14	9	5	5	15	8	12;
        2	4	7	9	5	4	9	11	2	14	5	14	11	16	16;
        2	4	13	12	12	7	15	10	5	2	15	5	13	7	4;
        2	5	9	9	3	12	8	14	15	12	14	5	3	2	15;
        2	5	11	7	2	11	9	4	16	7	16	9	14	14	4;
        2	6	2	13	3	3	12	9	7	16	6	9	16	13	12;
        2	6	9	7	7	16	13	3	12	2	13	12	9	16	6;
        2	7	12	15	2	12	4	10	13	15	13	4	5	5	10;
        2	7	14	16	5	9	2	9	16	11	11	5	7	4	14;
        2	8	5	12	5	2	14	14	8	15	3	9	12	15	9;
        2	9	13	4	2	13	8	11	6	4	6	8	15	15	11;
        2	10	3	2	13	16	8	10	8	13	11	11	16	3	5;
        2	11	15	3	11	6	14	10	15	10	6	7	7	14	3;
        2	16	4	5	16	14	7	11	4	11	14	9	9	7	5;
        3	3	4	6	11	12	13	6	12	14	4	5	13	5	14;
        3	3	6	5	16	9	15	5	9	10	6	4	15	4	10;
        3	4	5	14	4	6	12	13	5	13	6	11	11	12	14;
        3	4	9	16	10	4	16	15	3	5	10	5	15	6	6;
        3	4	16	10	5	10	4	9	9	16	15	6	3	5	15;
        3	5	12	11	14	5	11	13	3	6	14	6	13	4	4;
        3	6	4	10	6	5	9	15	4	15	5	16	16	9	10;
        3	7	8	8	16	11	12	4	15	11	4	7	16	3	15;
        3	7	16	11	4	15	3	15	11	12	12	4	7	8	16;
        3	8	7	15	4	8	15	12	3	16	4	16	12	11	11;
        3	8	15	4	16	4	8	7	7	15	12	11	3	16	12;
        3	10	10	15	16	5	4	6	16	4	3	15	9	6	9;
        3	13	11	5	4	12	4	11	6	6	5	3	14	13	12;
        3	14	7	9	14	10	13	8	7	8	10	4	4	13	9;
        5	5	8	14	16	13	6	14	13	7	8	15	6	15	7;
        5	6	11	7	10	8	5	8	7	12	12	10	6	9	11;
        5	6	13	8	13	5	7	7	6	16	14	15	8	16	15;
        5	7	9	10	7	11	6	12	9	12	11	8	8	6	10;
        5	9	6	8	10	9	8	12	5	11	10	11	12	7	7;
        5	10	10	12	8	11	9	7	8	9	5	12	6	7	6;
        5	10	12	6	5	12	8	9	7	6	7	8	11	11	9;
        5	13	15	15	14	8	6	7	16	8	7	13	14	5	16;
        9	10	13	10	11	15	15	9	16	12	14	13	16	14	11;
        9	11	12	15	12	9	13	13	11	14	10	16	15	14	16;
        9	12	10	15	13	14	9	14	15	11	11	13	12	16	10;];
    
%生成64*15相关矩阵；

for m=1:64
    for n=1:15
        for k=1:15
             switch   n
                case 1
                     correlation_value_64_15(m,n)=correlation_value_64_15(m,n)+correlation_value_16_15(table(m,k),k);
                 case 2
                     temp=table(m,1);
                    table(m,1:14)=table(m,2:15);
                    table(m,15)=temp;
                   correlation_value_64_15(m,n)=correlation_value_64_15(m,n)+correlation_value_16_15(table(m,k),k);
                case 3
                    temp=table(m,1);
                    table(m,1:14)=table(m,2:15);
                    table(m,15)=temp;
                    correlation_value_64_15(m,n)=correlation_value_64_15(m,n)+correlation_value_16_15(table(m,k),k);
                case 4
                    temp=table(m,1);
                    table(m,1:14)=table(m,2:15);
                    table(m,15)=temp;
                    correlation_value_64_15(m,n)=correlation_value_64_15(m,n)+correlation_value_16_15(table(m,k),k);
                case 5
                    temp=table(m,1);
                    table(m,1:14)=table(m,2:15);
                    table(m,15)=temp;
                    correlation_value_64_15(m,n)=correlation_value_64_15(m,n)+correlation_value_16_15(table(m,k),k);
                case 6
                    temp=table(m,1);
                    table(m,1:14)=table(m,2:15);
                    table(m,15)=temp;
                    correlation_value_64_15(m,n)=correlation_value_64_15(m,n)+correlation_value_16_15(table(m,k),k);
                case 7
                    temp=table(m,1);
                    table(m,1:14)=table(m,2:15);
                    table(m,15)=temp;
                    correlation_value_64_15(m,n)=correlation_value_64_15(m,n)+correlation_value_16_15(table(m,k),k);
                case 8
                    temp=table(m,1);
                    table(m,1:14)=table(m,2:15);
                    table(m,15)=temp;
                    correlation_value_64_15(m,n)=correlation_value_64_15(m,n)+correlation_value_16_15(table(m,k),k);
                case 9
                    temp=table(m,1);
                    table(m,1:14)=table(m,2:15);
                    table(m,15)=temp;
                    correlation_value_64_15(m,n)=correlation_value_64_15(m,n)+correlation_value_16_15(table(m,k),k);
                case 10
                    temp=table(m,1);
                    table(m,1:14)=table(m,2:15);
                    table(m,15)=temp;
                    correlation_value_64_15(m,n)=correlation_value_64_15(m,n)+correlation_value_16_15(table(m,k),k);
                case 11
                    temp=table(m,1);
                    table(m,1:14)=table(m,2:15);
                    table(m,15)=temp;
                    correlation_value_64_15(m,n)=correlation_value_64_15(m,n)+correlation_value_16_15(table(m,k),k);
                case 12
                   temp=table(m,1);
                    table(m,1:14)=table(m,2:15);
                    table(m,15)=temp;
                    correlation_value_64_15(m,n)=correlation_value_64_15(m,n)+correlation_value_16_15(table(m,k),k);
                case 13
                    temp=table(m,1);
                    table(m,1:14)=table(m,2:15);
                    table(m,15)=temp;
                    correlation_value_64_15(m,n)=correlation_value_64_15(m,n)+correlation_value_16_15(table(m,k),k);
                case 14
                    temp=table(m,1);
                    table(m,1:14)=table(m,2:15);
                    table(m,15)=temp;
                    correlation_value_64_15(m,n)=correlation_value_64_15(m,n)+correlation_value_16_15(table(m,k),k);
                case 15
                    temp=table(m,1);
                    table(m,1:14)=table(m,2:15);
                    table(m,15)=temp;
                    correlation_value_64_15(m,n)=correlation_value_64_15(m,n)+correlation_value_16_15(table(m,k),k);
             end
        end
    end
 end

%通过比较寻找出主扰码组号和帧起始位置；
A=0;
for m=1:64
    for n=1:15
        if (A<abs(correlation_value_64_15(m,n)))
             A=abs(correlation_value_64_15(m,n));
             r=m;
             c=n;
        end
    end
end
code_group_number =r-1
frame_boder=slot_boder+2560*(c-1);
phase_position=c;