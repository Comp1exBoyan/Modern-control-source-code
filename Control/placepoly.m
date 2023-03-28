function [K] = placepoly(inputA,inputB,dim,p,inputK)
syms s ;
fs = det(s*eye(dim)-inputA+inputB*inputK);
Q = coeffs(fs,s);
Q = flip(Q) ; 
f_s = 1;
for i = 1 : dim 
   f_s = f_s*(s-p(i)); 
end
P = coeffs(f_s,s);
P = flip(P);
for j = 1 : size(P)
eq(j) = P(j) - Q(j);
end
K = solve(eq,inputK);
disp(K);
end
