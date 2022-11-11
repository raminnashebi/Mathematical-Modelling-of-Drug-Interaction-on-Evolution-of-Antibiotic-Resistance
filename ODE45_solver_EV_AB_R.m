function dydt = Evolution_of_antibiotic_resistance(t,y,K,beta_s,beta_r,moo_s,moo_r,m,alpha11,alpha12,alpha21,alpha22,teta_1,teta_2,moo_1,moo_2,landa1,landa2)
% t : time
% y : vector which store bacterias population and antibiotic concentration 
% K      : Bacteria carrying capacity
% beta_s : The growth rate of sensitive bacteria
% beta_r : The growth rate of resistant bacteria
% moo_s: The natural death rate of sensitive bacteria
% moo_r: The natural death rate of resistant bacteria
% m :  The mutation rate of sensitive bacteria
% alpha11: effect of Antibiotic A on sensitive bacteria
% alpha12: effect of Antibiotic B on sensitive bacteria
% alpha21: effect of Antibiotic A on resistant bacteria
% alpha22: effect of Antibiotic B on resistant bacteria
% teta_1:
% teta_2:
% moo_1:
% moo_2:
% landa1:
% landa2:

 dydt = zeros(4,1);
 
  dydt(1) = beta_s*y(1)*(1-((y(1)+y(2))/K))-((m)*y(1))-(landa1*alpha11*alpha12*y(3)*y(4)+alpha11*y(3)+alpha12*y(4))*y(1)-(moo_s*y(1));
  dydt(2) = beta_r*y(2)*(1-((y(1)+y(2))/K))+((m)*y(1))-(landa2*alpha21*alpha22*y(3)*y(4)+alpha21*y(3)+alpha22*y(4))*y(2)-(moo_r*y(2));
  dydt(3)=teta_1-(moo_1*y(3));
  dydt(4)=teta_2-(moo_2*y(4));
   
 end
