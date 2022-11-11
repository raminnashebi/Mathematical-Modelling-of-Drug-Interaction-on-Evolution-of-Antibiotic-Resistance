
beta_s=1;
beta_r=0.65;
moo_s=0.5;
moo_r=0.5;
moo_1= (0.06/24)*50;
moo_2= (0.06/24)*50;
Emax_s=1.52*0.5;
Emax_r=1.12*0.5;
IC50_s=0.25*5;
IC50_r=5*5;



q1=(10^-6);
q2=(10^-8);
m=q1+q2;

gs=0;
landa1=1;
landa2=1;


alpha11=((-1 + sqrt(-(moo_s .* landa1)-(gs .* landa1)+(beta_s .* landa1)-(m.* landa1)+ 1)) ./ landa1);

teta_1=(alpha11*moo_1)/(Emax_s/IC50_s);
teta_2=(alpha11*moo_1)/(Emax_s/IC50_s);

%
alpha11=(Emax_s/IC50_s)*(moo_1/teta_1);
alpha12=(Emax_s/IC50_s)*(moo_2/teta_2);
alpha21=0.0367*alpha11;
alpha22=alpha21;


K=10^9;



tspan=[0 100];
y0=[0.9*10^9 0.1*10^9 0 0];
[t,y] = ode45(@(t,y)Evolution_of_antibiotic_resistance(t,y,K,beta_s,beta_r,moo_s,moo_r,m,alpha11,alpha12,alpha21,alpha22,teta_1,teta_2,moo_1,moo_2,landa1,landa2),tspan,y0);

RS=(beta_s/((landa1*alpha11*alpha12+alpha11+alpha12)+moo_s+q1+q2));
RR=(beta_r/((landa2*alpha21*alpha22+alpha21+alpha22)+moo_r));
m=q1+q2;
r=(m*((RS-1)/RS))/(beta_r*((1/RR)+(1/RS))+m);
s=((RS-1)/RS)-r;



ynew1(:,1)=(y(:,1)./(10^9));
ynew1(:,2)=(y(:,2)./10^9);
ynew1(:,3)=(y(:,3).*(moo_1/teta_1));
ynew1(:,4)=(y(:,4).*(moo_2/teta_2));


figure
semilogy(t,ynew1(:,1));
ylim([0.01 1])
hold on
semilogy(t,ynew1(:,2));
hold on
yyaxis right
semilogy(t,ynew1(:,3));
hold on
semilogy(t,ynew1(:,4));

