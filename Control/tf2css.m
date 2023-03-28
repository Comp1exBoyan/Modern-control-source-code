function [A,B,C] = tf2css(inputA,inputB,inputC,dim)
if dim <= 1
    disp("error!");
end
syms s ;
alpha = det(s * eye(dim) - inputA); 
K = sym2poly(alpha);
beta = zeros(dim);
A = zeros(dim);
B = zeros(dim,1);
C = zeros(1,dim);
for i = dim : -1 : 1
    beta(i) = inputC* inputA^(dim - i) * inputB;
    for j = 1 : dim - i
        beta(i) =  beta(i) + K(length(K)-(dim - j))*inputC*inputA^(dim - i-j)*inputB;
    end
end
    for j = 1 : dim 
       C(1,j)=beta(j); 
    end
for m = 1: dim-1
    for k = 1:dim
        if k == m +1
            A(m,k) = 1;
           
        end
        A(dim,k) = -K(length(K)-k + 1);
    end
end
B(dim,1) = 1;
end
