% Puntos equilibrio
% syms s r bs br alphas q ms mr mc u

syms s r

bs = 0.4;
br = 0.1;
alphas = 0.3960;
q = 4^(-3);
ms = 0.2;
mr = 0.09;
mc = 0.0083;
u = 0.12;

eq1 = (bs*s*( 1-(s+r) ) - alphas*s - ms*s + u);
eq2 = (br*r*( 1-(s+r) )+ q*alphas*s - mr*r );

peNum = vpasolve(eq1==0,eq2==0,s,r);
peSim = solve(eq1==0,eq2==0,s,r);

raices = roots([1 34699/84500 -54549/845000 29403/13520000]);

peSimEvS = [(1352000*raices(1)^2)/6039 + (632800*raices(1))/6039 - 480/61;
            (1352000*raices(2)^2)/6039 + (632800*raices(2))/6039 - 480/61;
            (1352000*raices(3)^2)/6039 + (632800*raices(3))/6039 - 480/61];
        
peSimEvR = [raices(1);
            raices(2);
            raices(3)];


% cond iniciales: [s, r, c] = [0.1,0,0]  MALO es [0,.1,0] (creo)


puntosEq = [peSimEvS peSimEvR];
% El punto de equilibrio es el segundo de puntosEq

% Puntos de equilibrio hallado experimentalmente
% so = 0.3356; ro = 0.0683 ; uo= 0.12
so = 0.3356; ro = 0.0683 ; uo= 0.12;
%Cambio punto en el cual esta linealizado y sistema utilizado
[A,B,C,D] = linmod('resistenciaBacterias_entrada',[so,ro],uo);
eig(A);
% so=puntosEq(1,1);ro=puntosEq(1,2);co=puntosEq(1,3);

% A = [-0.4918,-0.1342,-0.0006,-0.0372]
%Matrices modelo lineal teoricamente (JACOBIAN)
At = [bs*(1-2*so-ro) - alphas - ms, - bs*so ;
     -br*ro + q*alphas, br*(1-2*ro-so) - mr ];
Bt = [1 0 ]';
Ct = [1 1 ; 0 1 ];
Dt = [0 ]';


%Comparacion mediante simulink de el modelo lineal vs el no lineal 
% Metodo grafico
tiempo = tout;
plot(tiempo,SRout,'r-',tiempo,Rout,'b-',tiempo,SRLout,'r:',tiempo,RLout,'b:')
title('Comparacion modelo lineal con modelo no lineal u = 0.01')
ylabel('Adimensional')
xlabel('Tiempo(min)')
legend('y_1-no lineal','y_2-no lineal','y_1-(Lineal)','y_2-(Lineal)')

 matlab2tikz('comparacion_u04.tex')