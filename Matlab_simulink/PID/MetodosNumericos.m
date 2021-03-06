%metodo euler

ts = 5; % delta t
tiempoFinal = 700;

% Si quisiera hacer el mismo numero de puntos que la simulacion
% n = size(tout);
% n = n(1);

% Numero de puntos fijo
n = ceil(tiempoFinal/ts);

se = zeros(1,n);
re = zeros(1,n);
ce = zeros(1,n);
% vector tiempo
tiempo = (0:ts:ts*(n-1));

% Condiciones iniciales
se(1) = 0.2; re(1) = 0; ce(1) = 0;

% Parametros
param = [0.14 0.1 0.3960 4^(-3) 0.2 0.09 0.0083];
param = num2cell(param);
[bs,br,alphas,q,ms,mr,mc] = param{:};

% ciclo metodo de euler
for k = 1:(n-1)
    
    s = se(k);
    r = re(k);
    c = ce(k);
    se(k+1) = s + (bs*s*( 1-(s+r) ) - alphas*c*s - ms*s)*ts;
    re(k+1) = r + (br*r*( 1-(s+r) )+ q*alphas*c*s - mr*r )*ts;
    ce(k+1) = c + (mc - mc*c)*ts;
end

%Runge-Kutta

sr = zeros(1,n);
rr = zeros(1,n);
cr = zeros(1,n);

sr(1) = 0.2; rr(1) = 0; cr(1) = 0;

ds = @(s,r,c) (bs*s*( 1-(s+r) ) - alphas*c*s - ms*s);
dr = @(s,r,c) (br*r*( 1-(s+r) )+ q*alphas*c*s - mr*r );
dc = @(s,r,c) (mc - mc*c);

for i=1:(n-1)
   s = sr(i);
   r = rr(i);
   c = cr(i);
    
   k1 = ds(s,r,c);
   m1 = dr(s,r,c);
   n1 = dc(s,r,c);
   
   k2 = ds(s + ts*k1/2, r+ ts*m1/2, c+ ts*n1/2);  
   m2 = dr(s + ts*k1/2, r+ ts*m1/2, c+ ts*n1/2); 
   n2 = dc(s + ts*k1/2, r+ ts*m1/2, c+ ts*n1/2); 
   
   k3 = ds(s + ts*k2/2, r+ ts*m2/2, c+ ts*n2/2);  
   m3 = dr(s + ts*k2/2, r+ ts*m2/2, c+ ts*n2/2); 
   n3 = dc(s + ts*k2/2, r+ ts*m2/2, c+ ts*n2/2); 
   
   k4 = ds(s+k3*ts, r+m3*ts, c+n3*ts);
   m4 = dr(s+k3*ts, r+m3*ts, c+n3*ts);
   n4 = dc(s+k3*ts, r+m3*ts, c+n3*ts);
   
   sr(i+1) = s + (k1+(k2*2)+(k3*2)+k4)*ts/6;
   rr(i+1) = r+ (m1+(m2*2)+(m3*2)+m4)*ts/6;
   cr(i+1) = c + (n1+(n2*2)+(n3*2)+n4)*ts/6;
end

%Symulink
load_system('resistenciaBacterias');
set_param('resistenciaBacterias/Subsystem','par','[0.14 0.1 0.3960 4^(-3) 0.2 0.09 0.0083]','x0','[0.2 0 0]');
set_param('resistenciaBacterias','StopTime',num2str(tiempoFinal));
sim('resistenciaBacterias');


plot(tiempo,re,'r:',tiempo,se+re,'b:',tiempo,rr,'r-',...
   tiempo,sr+rr,'b-',tout,Rout,'r--',tout,Sout + Rout,'b--')
title('Comparación salidas métodos númericos')
ylabel('Adimensional')
xlabel('Tiempo(min)')
legend('y1-(euler)','y2-(euler)',...
    'y1-(Runge-Kutta)','y2-(Runge-Kutta)'...
    ,'y1-(Symulink)','y2-(Symulink)')

% plot(tiempo,se,'r:',tiempo,re,'b:',tiempo,ce,'k:',tiempo,sr,'r-',...
%    tiempo,rr,'b-',tiempo,cr,'k-',tout,Sout,'r--',tout,Rout,'b--',tout,Cout,'k--')
% % plot(tiempo,se,'r:',tiempo,re,'r-',tiempo,ce,'r--',tiempo,sr,'b:',...
% %     tiempo,rr,'b-',tiempo,cr,'b--',tout,Sout,'k:',tout,Rout,'k-',tout,Cout,'k--')
% 
% ylabel('Adimensional')
% xlabel('Tiempo(min)')
% legend('Tasa bact. sensibles-(euler)','Tasa bact. resistentes-(euler)','Concent. antibiótico-(euler)',...
%     'Tasa bact. sensibles-(Runge-Kutta)','Tasa bact. resistentes-(Runge-Kutta)','Concent. antibiótico-(Runge-Kutta)'...
%     ,'Tasa bact. sensibles-(Symulink)','Tasa bact. resistentes-(Symulink)','Concent. antibiótico-(Symulink)')