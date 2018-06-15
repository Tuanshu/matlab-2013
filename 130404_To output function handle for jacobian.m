clear all

syms Q1 Q2 Q3 Q4 Q5 Q6 Q7 Q8 Q9 Lambda

P=[Q1 Q2 Q3; Q4 Q5 Q6; Q7 Q8 Q9];
I=[1 0 0; 0 1 0; 0 0 1];
PT=P.';
Pi=inv(P);


P_LA=simplify(((P+Lambda*(PT\I))));
P_LMA=simplify(((P+Lambda*(PT\[P(1,1) 0 0; 0 P(2,2) 0; 0 0 P(2,2)]))));