
beta_s=1;
beta_r=0.65;
moo_s=0.5;
moo_r=0.5;
moo_1= (0.06/24)*50;%(0.06/24)*50;
moo_2= (0.06/24)*50;%(0.06/24)*50;%(0.05/24);%*84;

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
%{
teta_1=0.060216386643907
teta_2=0.046205686484528
%}

%{
teta_1=1.204327732878142e-04;
teta_2=9.241137296905680e-05;
teta_3=1.204327732878142e-04;
teta_4=9.241137296905680e-05;
%}



%{
alpha11=((-1 + sqrt(-(moo_s .* landa1)-(gs .* landa1)+(beta_s .* landa1)-(m.* landa1)+ 1)) ./ landa1);
alpha12=alpha11;
alpha21=0.0367*alpha11;
alpha22=alpha21;
%}

%
alpha11=(Emax_s/IC50_s);%*(moo_1/teta_1);
alpha12=(Emax_s/IC50_s);%*(moo_2/teta_2);
alpha21=0.0367*alpha11;
alpha22=alpha21;
%alpha21=(Emax_r/IC50_r);%*(moo_1/teta_1);
%alpha22=(Emax_r/IC50_r);%*(moo_2/teta_2);
%}

K=10^9;



tspan=[0 100];
%y0=[6000 20 4 4];
y0=[0.9*10^9 0.1*10^9 0 0];
[t,y] = ode45(@(t,y)Evolution_of_antibiotic_resistance(t,y,K,beta_s,beta_r,moo_s,moo_r,m,alpha11,alpha12,alpha21,alpha22,teta_1,teta_2,moo_1,moo_2,landa1,landa2),tspan,y0);

RS=(beta_s/((landa1*alpha11*alpha12+alpha11+alpha12)+moo_s+q1+q2));
RR=(beta_r/((landa2*alpha21*alpha22+alpha21+alpha22)+moo_r));
m=q1+q2;
r=(m*((RS-1)/RS))/(beta_r*((1/RR)+(1/RS))+m);
s=((RS-1)/RS)-r;

%{
%figure
plot(t,y(:,1));
hold on
plot(t,y(:,2));
%xlim([0 10])
%ylim([0 1200])
hold on
%yyaxis right
%plot(t,y(:,3));
hold on
%plot(t,y(:,4));
%yyaxis left
%}

ynew1(:,1)=(y(:,1)./(10^9));
ynew1(:,2)=(y(:,2)./10^9);
ynew1(:,3)=(y(:,3).*(moo_1/teta_1));
ynew1(:,4)=(y(:,4).*(moo_2/teta_2));


figure
%semilogy(t,ynew1(:,1));
%ylim([0.01 1])
hold on
semilogy(t,ynew1(:,2));
hold on
yyaxis right
semilogy(t,ynew1(:,3));
hold on
semilogy(t,ynew1(:,4));




x=[-1.667 -1.65 -1.6 -1.517 -1.4 -1.25 -1.067 -0.85 -0.6 -0.3167 0 0.35 0.7333 1.15 1.6];
G1 = 1:0.1:2.4;               

plot(x, G1)
hold on 
patch([x fliplr(x)], [G1 min(ylim)*ones(size(G1))], 'g')
patch([x fliplr(x)], [G1 max(ylim)*ones(size(G1))], 'r')
                                                              % Below Lower Curve
hold off










[X Z]=meshgrid(y(:,3),y(:,4));

Rs1=(1-(1./(1+(((alpha11.*X))+((alpha12.*Z))+((landa1).*((alpha11*alpha12.*X.*Z)))))));
Rs2=(1-(1./(1+(((alpha21.*X))+((alpha22.*Z))+((landa2).*((alpha21*alpha22.*X.*Z)))))));


subplot(1,3,2)
surf(X,Z,Rs1)

subplot(1,3,3)
surf(X,Z,Rs2)


figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create surf
surf(X,Z,Rs1,'Parent',axes1,'EdgeColor','none');

view(axes1,[0.599999999999994 90]);
grid(axes1,'on');


figure2 = figure;

% Create axes
axes2 = axes('Parent',figure2);
hold(axes2,'on');

% Create surf
surf(X,Z,Rs2,'Parent',axes2,'EdgeColor','none');

view(axes2,[0.599999999999994 90]);
grid(axes2,'on');
