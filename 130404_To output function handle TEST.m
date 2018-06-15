clear all


t1= @ (n) (2*n1_Considered)./(n1_Considered+n); 
t1_r= @ (n) (2*n)./(n2_Considered+n);
t2= @ (n) (2*n)./(n+n2_Considered); 
r1= @ (n) (n1_Considered-n)./(n1_Considered+n);
r1_r= @ (n) (n-n1_Considered)./(n+n1_Considered);
r2= @ (n) (n-n2_Considered)./(n+n2_Considered);

A = @ (n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now) abs(Ratio_Upper2Reference_Now./r_AN100_Considered.*(r1(n)));
B = @ (n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now) abs(Ratio_Upper2Reference_Now./r_AN100_Considered.*(Ratio_Lower2Upper_Now.*(t1(n).*t1_r(n).*r2(n).*exp(1i*4*pi.*Frequency_Considered.*n.*Thickness_Now/c))));
C = @ (n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now) diff(unwrap(angle(1./r1(n).*(Ratio_Lower2Upper_Now.*(t1(n).*t1_r(n).*r2(n).*exp(1i*4*pi.*Frequency_Considered.*n.*Thickness_Now/c)))),[],1)); %雖然沒有絕對的phase, 但是是不是應該只shift pi的整數?

Q1_P_add_Lambda_PT_inv = @ (Q1,Q2,Q3,Q4,Q5,Q6,q7,Q8,Q9,Lambda) Q1 + (Q5.*Q9 - Q6.*Q8)./(Lambda*(Q1.*Q5.*Q9 - Q1.*Q6.*Q8 - Q2.*Q4.*Q9 + Q2.*Q6.*Q7 + Q3.*Q4.*Q8 - Q3.*Q5.*Q7));
Q2_P_add_Lambda_PT_inv = @ (Q1,Q2,Q3,Q4,Q5,Q6,q7,Q8,Q9,Lambda) Q2 - (Q4.*Q9 - Q6.*Q7)./(Lambda*(Q1.*Q5.*Q9 - Q1.*Q6.*Q8 - Q2.*Q4.*Q9 + Q2.*Q6.*Q7 + Q3.*Q4.*Q8 - Q3.*Q5.*Q7));
Q3_P_add_Lambda_PT_inv = @ (Q1,Q2,Q3,Q4,Q5,Q6,q7,Q8,Q9,Lambda) Q3 + (Q4.*Q8 - Q5.*Q7)./(Lambda*(Q1.*Q5.*Q9 - Q1.*Q6.*Q8 - Q2.*Q4.*Q9 + Q2.*Q6.*Q7 + Q3.*Q4.*Q8 - Q3.*Q5.*Q7));
Q4_P_add_Lambda_PT_inv = @ (Q1,Q2,Q3,Q4,Q5,Q6,q7,Q8,Q9,Lambda) Q4 - (Q2.*Q9 - Q3.*Q8)./(Lambda*(Q1.*Q5.*Q9 - Q1.*Q6.*Q8 - Q2.*Q4.*Q9 + Q2.*Q6.*Q7 + Q3.*Q4.*Q8 - Q3.*Q5.*Q7));
Q5_P_add_Lambda_PT_inv = @ (Q1,Q2,Q3,Q4,Q5,Q6,q7,Q8,Q9,Lambda) Q5 + (Q1.*Q9 - Q3.*Q7)./(Lambda*(Q1.*Q5.*Q9 - Q1.*Q6.*Q8 - Q2.*Q4.*Q9 + Q2.*Q6.*Q7 + Q3.*Q4.*Q8 - Q3.*Q5.*Q7));
Q6_P_add_Lambda_PT_inv = @ (Q1,Q2,Q3,Q4,Q5,Q6,q7,Q8,Q9,Lambda) Q6 - (Q1.*Q8 - Q2.*Q7)./(Lambda*(Q1.*Q5.*Q9 - Q1.*Q6.*Q8 - Q2.*Q4.*Q9 + Q2.*Q6.*Q7 + Q3.*Q4.*Q8 - Q3.*Q5.*Q7));
Q7_P_add_Lambda_PT_inv = @ (Q1,Q2,Q3,Q4,Q5,Q6,q7,Q8,Q9,Lambda) Q7 + (Q2.*Q6 - Q3.*Q5)./(Lambda*(Q1.*Q5.*Q9 - Q1.*Q6.*Q8 - Q2.*Q4.*Q9 + Q2.*Q6.*Q7 + Q3.*Q4.*Q8 - Q3.*Q5.*Q7));
Q8_P_add_Lambda_PT_inv = @ (Q1,Q2,Q3,Q4,Q5,Q6,q7,Q8,Q9,Lambda) Q8 - (Q1.*Q6 - Q3.*Q4)./(Lambda*(Q1.*Q5.*Q9 - Q1.*Q6.*Q8 - Q2.*Q4.*Q9 + Q2.*Q6.*Q7 + Q3.*Q4.*Q8 - Q3.*Q5.*Q7));
Q9_P_add_Lambda_PT_inv = @ (Q1,Q2,Q3,Q4,Q5,Q6,q7,Q8,Q9,Lambda) Q9 + (Q1.*Q5 - Q2.*Q4)./(Lambda*(Q1.*Q5.*Q9 - Q1.*Q6.*Q8 - Q2.*Q4.*Q9 + Q2.*Q6.*Q7 + Q3.*Q4.*Q8 - Q3.*Q5.*Q7));

det_array = @ (Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9) Q1.*Q5.*Q9 + Q2.*Q6.*Q7 + Q3.*Q4.*Q8 -Q3.*Q5.*Q7 -Q2.*Q4.*Q9 -Q1.*Q6.*Q8;
 
T_abs = @ (n,Thickness_Now) abs(t_AN100_Considered.*t1(n).*t2(n).*exp(1i*2*pi*Frequency_Considered.*n.*Thickness_Now/c)./(1-r1_r(n).*r2(n).*exp(1i*4*pi*Frequency_Considered.*n.*Thickness_Now/c))).^2;          %Note, 120827 remove 2 from 分子4*pi > 2*pi
                    

syms n n1_Considered n2_Considered Thickness_Now Ratio_Lower2Upper_Now Ratio_Upper2Reference_Now Frequency_Considered r_AN100_Considered


t1= (2*n1_Considered)./(n1_Considered+n); 
t1_r=  (2*n)./(n2_Considered+n);
t2=  (2*n)./(n+n2_Considered); 
r1=  (n1_Considered-n)./(n1_Considered+n);
r1_r=  (n-n1_Considered)./(n+n1_Considered);
r2=  (n-n2_Considered)./(n+n2_Considered);

A=abs(Ratio_Upper2Reference_Now./r_AN100_Considered.*(r1));
B = abs(Ratio_Upper2Reference_Now./r_AN100_Considered.*(Ratio_Lower2Upper_Now.*(t1.*t1_r.*r2.*exp(1i*4*pi.*Frequency_Considered.*n.*Thickness_Now/c))));
C = diff(unwrap(angle(1./r1.*(Ratio_Lower2Upper_Now.*(t1.*t1_r.*r2.*exp(1i*4*pi.*Frequency_Considered.*n.*Thickness_Now/c)))),[],1)); 

JACA=jacobian(P,X);

JACA=matlabFunction(JACA);








Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9
dn_real = det_array(delta_A,Q2,Q3,delta_B,Q5,Q6,delta_C,Q8,Q9)./det_array(Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9);
dn_imag = det_array(Q1,delta_A,Q3,Q4,delta_B,Q6,Q7,delta_C,Q9)./det_array(Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9);
dthickness = det_array(Q1,Q2,delta_A,Q4,Q5,delta_B,Q7,Q8,delta_C)./det_array(Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9);
                        