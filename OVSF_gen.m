function y = OVSF_gen(x)
[m,n] = size(x);
y = zeros(2*n);
for i = 1:n
    y(2*i-1:2*i,:) = [x(i,:),x(i,:);x(i,:),-x(i,:)];
end