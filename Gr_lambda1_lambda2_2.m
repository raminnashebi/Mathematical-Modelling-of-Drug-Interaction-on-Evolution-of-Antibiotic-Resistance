beta_s=1;
beta_r=0.65;
moo_s=0.5;
moo_r=0.5;
teta_1=0.21;
teta_2=0.21;
moo_1= 0.012;
moo_2= 0.012;
Emax_s=1.52;
Emax_r=1.12;
IC50_s=0.25;
IC50_r=5;



q1=(10^-6);
q2=(10^-8);
m=q1+q2;


alpha11=(Emax_s/IC50_s)*(moo_1/teta_1);
alpha12=(Emax_s/IC50_s)*(moo_2/teta_2);
alpha21=(Emax_r/IC50_r)*(moo_1/teta_1);
alpha22=(Emax_r/IC50_r)*(moo_2/teta_2);

K=0.0367;
gs=0;

syms lambda1 lambda2
lambda1=-1.66:0.1:1.66;
lambda2=-1.66:0.1:1.66;
[lambda1 lambda2]=meshgrid(lambda1,lambda2);
gr = beta_r - ((K^2 .* ((-1 + sqrt(-(moo_s .* lambda1)-(gs .* lambda1) + (beta_s .* lambda1)-(m.* lambda1) + 1)) .^ 2) .* lambda2 )./ lambda1 .^ 2) - (((2*K) .* (-1 + sqrt(-(moo_s .* lambda1)-(gs .* lambda1)+(beta_s .* lambda1)-(m.* lambda1)+ 1))) ./ lambda1)-moo_r;
figure
surf(lambda1,lambda2,gr)
