function [u_1, u_2] = MatrixShrink(A, B, c, d_1, d_2, alpha) 
% "MATRIXSHRINK"
% Shrinkage as in bachelor thesis, compoutes argmin \| A u_1 + B u_2 - c \|_1 + alpha/2 \| u - d \|_2^2
% written by:       Jan Henrik Fitschen         11/03/2014

case_one 	= (~A)&(~B);
case_two 	= (A)&(~B);
case_three 	= (~A)&(B);
%case_four 	= (A)&(B);
c_A = c ./ (A+eps);
c_B = c ./ (B+eps);

u_1_two   = SoftShrink(-c_A + d_1, 1/alpha * abs(A)) + c_A;
u_2_three = SoftShrink(-c_B + d_2, 1/alpha * abs(B)) + c_B;
            
d1 = A .* d_1;
d2 = B .* d_2 - c;
lambda1 = 1/alpha * A .* A;
lambda2 = 1/alpha * B .* B;
            
t1 = ((d1+d2-lambda1-lambda2) > 0);
t2 = ((d1+d2+lambda1+lambda2) < 0);
t3 = ((~t1) & (~t2));
x = (d1 - lambda1) .* t1 + ...
    (d1 + lambda1) .* t2 + ...
    (d1 .* lambda2 - d2 .* lambda1)./(lambda1+lambda2) .* t3;
y = (d2 - lambda2) .* t1 + ...
    (d2 + lambda2) .* t2 + ...
    (d2 .* lambda1 - d1 .* lambda2)./(lambda1+lambda2) .* t3;
u_1_four = x ./ A;
u_2_four = (y + c) ./ B;

u_1             = u_1_four;             % case four is the most like case, this is faster than
u_2             = u_2_four;             % u_1(case_four) = u_1_four(case_four); 
u_1(case_one) 	= d_1(case_one);
u_2(case_one) 	= d_2(case_one);
u_1(case_two) 	= u_1_two(case_two);
u_2(case_two) 	= d_2(case_two);
u_1(case_three) = d_1(case_three);
u_2(case_three) = u_2_three(case_three);

end
