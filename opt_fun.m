function y = opt_fun(x)

global low_index upp_index gamma_xy_fit w L_fit Phi_fit skip

%s = tf('s');

Wg = 1; Wp = 1;
       
%West = tf(x(1),[1/x(2) 1]);
West = pretpostavljena(x);

[AMP,PHI] = bode(West,w(low_index:skip:upp_index));
L = 20*log10(AMP(1,:)');
PHI = PHI(1,:)'*pi/180;
ll = length(low_index:skip:upp_index);

a =  sum( ( (1.58  * (1 - exp( - (gamma_xy_fit).^2))).^2 .*  (   L_fit - L    ).^2 ) )*Wg; 
b =  (180/pi)*sum( ( (1.58  * (1 - exp( - (gamma_xy_fit).^2))).^2 .*  ( Phi_fit - PHI   ).^2 ) )*Wp;  

y = 20/(ll)*( a + b );  


end