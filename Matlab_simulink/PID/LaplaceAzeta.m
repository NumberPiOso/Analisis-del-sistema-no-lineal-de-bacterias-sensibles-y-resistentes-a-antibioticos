%(vpa(inv(s*eye(2) - A),5))
%vpa((C*inv(s*eye(2) - A)*B + D),5)
%tf(sys)
%Laplace a transf Z

so = 0.3356; ro = 0.0683 ; uo= 0.12;
%Cambio punto en el cual esta linealizado y sistema utilizado
[A,B,C,D] = linmod('resistenciaBacterias_entrada',[so,ro],uo);

% Transformada Z
Ts = 0.2;
% t = k*Ts;

%Sistema lineal
sys = ss(A,B,C,D);
TrasG = tf(sys);
TrasZ = tf(sys,Ts);
%Funcion de transferencia (hallada en otro codigo)
syms s z k
%G1 = (s + 0.0062)/(s^2 + .2460*s + .00236);
%G2 = (-0.0038)/(s^2 + .2460*s + .00236);

%roots([1 .2460 .00236]) raices polinomios

% Discretizacion de transformada laplace a z 
G1 = (s + 0.03658)/ (s^2 + 0.529 *s + 0.01822);
G2 = -0.0006425/ (s^2 + 0.529 *s + 0.01822);
% Hallo la inversa G(S)/s y cambio t por K*Ts
Art = eval(ilaplace(G1/s)); 
GZ1 = (1-z^(-1))*ztrans(Art)
GZ2 = (1-z^(-1))*ztrans(eval(ilaplace(G2/s)))
% Discretizacion por matlab (para verificar)
sysD = c2d(sys,Ts);
tdd =  tf(sysD);


% Ahora si empieza lo bueno bueno, comparaciones

%Comparacion salida SR
% SRLout=SRLout(1:length(SRztrans));
% tout =tout(1:length(SRztrans));
% SRltrans = SRltrans(1:length(SRztrans));
% plot(tout ,SRLout ,tout ,SRztrans,tout,SRltrans)
% title('Comparación y1 transformadas vs lineal')
% ylabel('Adimensional')
% xlabel('Tiempo(min)')
% legend('y1-(lineal)','y1-(transf. Z)','y1-(transf. Laplace)')

%COmparacion 
% RLout=RLout(1:length(Rztrans));
% tout =tout(1:length(Rztrans));
% Rltrans = Rltrans(1:length(Rztrans));
% plot(tout ,RLout ,tout ,Rztrans,tout,Rltrans)
% title('Comparación y2 transformadas vs lineal')
% ylabel('Adimensional')
% xlabel('Tiempo(min)')
% legend('y2-(lineal)','y2-(transf. Z)','y2-(transf. Laplace)')


%codigos utilizados para latex en orden
%latex(vpa(partfrac(G1/s),5))
% latex(vpa(ilaplace(G1/s),5))
%latex(vpa(eval(ilaplace(G1/s)),5))
% vpa(simplifyFraction(GZ1),5) y lo paso manualmente, (bien gay) a decimal
            % a1 = .00074996/(1.0001*1.0024);
            % a1*13.351
            % a1*13.35
%manual = latex((0.01*z +0.01)/((z-.9999)*(z-.9976)))
% numerador G2(z) 1.8984e-07z  - 1.8984e-07
a2 = -0.014249/(1.0001*1.0024);
a2*0.000013356;
a2*0.000013345;

% ============================================================================
% Reducción de Polos Transformadas continuas
G1r = 0.9878*1/(s + 0.4920);
G2r = 2.0408 *(-0.00064)/(s + 0.0370);


% Entrada Seno para el sistema
tsin = 0:Ts:900;  
usin = sin(tsin);  
% =============================================================================
% Reducción por cancelación
% Con tf de la 1
sys1 = tf([1 0.0365],[1 0.529 0.01822]);
sysr1 = tf(0.9878,[1 0.4920]);
       % Escalon 
figure(1)
hold on; step(sys1); step(sysr1);legend('Sistema original TF1','Sistema reducido TF1') ;hold off
       % Seno                                  
sys1Sin = lsim(sys1, usin, tsin);     
sysr1Sin = lsim(sysr1, usin, tsin);

figure(2)
hold on; plot(tsin,sys1Sin);plot(tsin,sysr1Sin);
title('Sin Response')
ylabel('Amplitud')
xlabel('Tiempo')
legend('Sistema original TF1','Sistema reducido TF1');hold off;

% Con tf de la 2
sys2 = tf(-0.00064,[1 0.529 0.01822]);
sysr2 = tf(2.0408 *(-0.00064),[1 0.0370]);
    % Escalon
figure(3)
hold on; step(sys2); step(sysr2);legend('Sistema original TF2','Sistema reducido TF2') ;hold off
    % Seno
sys2Sin = lsim(sys2, usin, tsin);     
sysr2Sin = lsim(sysr2, usin, tsin);

figure(4)
hold on; plot(tsin,sys2Sin);plot(tsin,sysr2Sin);
title('Sin Response')
ylabel('Amplitud')
xlabel('Tiempo')
legend('Sistema original TF2','Sistema reducido TF2');hold off;
                          
% ==========================================================================

% Reducción por reduce

sys1Reduce = reduce(sys1,1);
sys2Reduce = reduce(sys2,1);

%Tf1
    % Escalon
figure(5)
hold on; step(sys1); step(sys1Reduce);legend('Sistema original TF1','Sistema reducido TF1 (reduce)') ;hold off

% Seno                                  
sys1Sinreduce = lsim(sys1, usin, tsin);     
sysr1Sinreduce = lsim(sys1Reduce, usin, tsin);

figure(6)
hold on; plot(tsin,sys1Sinreduce);plot(tsin,sysr1Sinreduce);
title('Sin Response')
ylabel('Amplitud')
xlabel('Tiempo')
legend('Sistema original TF1','Sistema reducido TF1');hold off;

% Tf2
figure(7)
hold on; step(sys2); step(sys2Reduce);legend('Sistema original TF2','Sistema reducido TF2 (reduce)') ;hold off
% Seno                                  
sys2Sinreduce = lsim(sys2, usin, tsin);     
sysr2Sinreduce = lsim(sys2Reduce, usin, tsin);

figure(8)
hold on; plot(tsin,sys2Sinreduce);plot(tsin,sysr2Sinreduce);
title('Sin Response')
ylabel('Amplitud')
xlabel('Tiempo')
legend('Sistema original TF1','Sistema reducido TF1');hold off;

% =============================================================================
% Balred

sys1balred = balred(sys1,1);
sys2balred = balred(sys2,1);


%Tf1
    % Escalon
figure(9)
hold on; step(sys1); step(sys1balred);legend('Sistema original TF1','Sistema reducido TF1 (balred)') ;hold off

% Seno                                  
sys1Sinbalred = lsim(sys1, usin, tsin);     
sysr1Sinbalred = lsim(sys1balred, usin, tsin);

figure(10)
hold on; plot(tsin,sys1Sinbalred);plot(tsin,sysr1Sinbalred);
title('Sin Response')
ylabel('Amplitud')
xlabel('Tiempo')
legend('Sistema original TF1','Sistema reducido TF1');hold off;

% Tf2
figure(11)
hold on; step(sys2); step(sys2balred);legend('Sistema original TF2','Sistema reducido TF2 (balred)') ;hold off
% Seno                                  
sys2Sinbalred = lsim(sys2, usin, tsin);     
sysr2Sinbalred = lsim(sys2balred, usin, tsin);

figure(12)
hold on; plot(tsin,sys2Sinbalred);plot(tsin,sysr2Sinbalred);
title('Sin Response')
ylabel('Amplitud')
xlabel('Tiempo')
legend('Sistema original TF1','Sistema reducido TF1');hold off;


% =============================================================================

% minreal

sys1minreal = minreal(sys1,1);
sys2minreal = minreal(sys2,1);


%Tf1
    % Escalon
figure(13)
hold on; step(sys1); step(sys1minreal);legend('Sistema original TF1','Sistema reducido TF1 (minreal)') ;hold off

% Seno                                  
sys1Sinminreal = lsim(sys1, usin, tsin);     
sysr1Sinminreal = lsim(sys1minreal, usin, tsin);

figure(14)
hold on; plot(tsin,sys1Sinminreal);plot(tsin,sysr1Sinminreal);
title('Sin Response')
ylabel('Amplitud')
xlabel('Tiempo')
legend('Sistema original TF1','Sistema reducido TF1');hold off;

% Tf2
figure(15)
hold on; step(sys2); step(sys2minreal);legend('Sistema original TF2','Sistema reducido TF2 (minreal)') ;hold off
% Seno                                  
sys2Sinminreal = lsim(sys2, usin, tsin);     
sysr2Sinminreal = lsim(sys2minreal, usin, tsin);

figure(16)
hold on; plot(tsin,sys2Sinminreal);plot(tsin,sysr2Sinminreal);
title('Sin Response')
ylabel('Amplitud')
xlabel('Tiempo')
legend('Sistema original TF1','Sistema reducido TF1');hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cuarto Punto
syms z
% Funciones de transferencia discretas del modelo
GZ1 = (0.1905*z - 0.1891)/( z^2 - 1.899*z + 0.8996);
GZ2 = (-1.241e-05*z - 1.198e-05)/(z^2 - 1.899*z + 0.8996);

% secuencia de ponderacion para ft1 y ft2, 10 valores
sp1 = zeros(10,1);
sp2 = zeros(10,1);

for i = 1:10
    % calcula la transformada inversa z y la evalua en i - 1 para que
    % empiece en 0
    sp1(i) =  iztrans(GZ1,(i-1));
    sp2(i) =  iztrans(GZ2,(i-1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



