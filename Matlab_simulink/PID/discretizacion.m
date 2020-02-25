% Puntos de linealizacion
uo = 0.12;
so = 0.2;
ro = 0;
tiempoFinal = 800;
n = 25;
yss = zeros(4,n);
%Symulink
load_system('VariablesEstado');
% set_param('VariablesEstado/Subsystem','par','[0.14 0.1 0.3960 4^(-3) 0.2 0.09 0.0083]','x0','[0.2 0 0]');
set_param('VariablesEstado','StopTime',num2str(tiempoFinal));
for k = 1:n
    so = so + 0.1;
    ro = ro + 0.1;
    sim('VariablesEstado');
    yss(1,k) = mean(Sout(length(Sout)-20:length(Sout)));
    yss(2,k) = mean(Rout(length(Rout)-20:length(Rout)));
    yss(3,k) = so;
    yss(4,k) = ro;
end

plot(tsim,Rout,tsim,Sout)


%Discretizacion
Ts = 0.5;

so = 0.335599308431370; ro = 0.068323868592545; uo= 0.12;
[A,B,C,D] = linmod('resistenciaBacteriasSalRes_entrada',[so,ro],uo);
sys = ss(A,B,C,D);
sysD = c2d(sys,Ts);


%  step(sys(1))
%
%sys = sim('resistB_Compara','SaveState','on','StateSaveName','xout','SaveOutput','on','OutputSaveName','yout');

%--------------------------------------------------------
% Asignacion de polas
saturMax = 2;
saturMin = -.5;
%[B,A] = numden(tf(sysD));
syms qs0 qs1 qs2 rs1 z 
Bp = (-0.003178*z - 0.002949);
Ap = (z^2 -1.794*z + 0.799);
Qc = qs0*z^2 + qs1*z + qs2;
Pc = (z-rs1)*(z-1);

%PolosDes = [.98 .99 -.999 -.95];
PolosDes = [0 0 0 0];
Ad = (z-PolosDes(1))*(z-PolosDes(2))*(z-PolosDes(3))*(z-PolosDes(4));
vecAsig= coeffs(Ap*Pc + Bp*Qc - Ad,z);
solucionParam = solve(vecAsig==[0 0 0 0],qs0,qs1,qs2,rs1)

q0 = eval(solucionParam.qs0);q1 = eval(solucionParam.qs1);
q2 = eval(solucionParam.qs2);r1 = eval(solucionParam.rs1);
%-----------------------------------------------------
% Matriz de controlabilidad

Mc = [B A*B];
cond(Mc)

Co = ctrb(sysD.A,sysD.B);
cond(Co)
%-----------------------------------------------------
% Realimentación 


get(0,'Factory')
set(0,'defaultfigurecolor',[1 1 1])
set(groot,'defaultLineLineWidth',4)
set(groot,'DefaultAxesFontSize', 20)
%set(groot,'Backgrouncolor', 'w')

p = [0 0];
%K = place(sysD.A,sysD.B,[p(1) p(2)]);
K = acker(sysD.A,sysD.B,[p(1) p(2)]);
sim('RealNoLineal')



plot(tsim,Rout)
title('Baterias resistentes')
xlabel('tiempo (minutos)')
ylabel('Bacterias resistentes (adim.)')


plot(tsim,ControlAccion)
title('Acción de control')
xlabel('tiempo (minutos)')
ylabel('Acción control (Adim)')

plot(tsim,perturbacion)
title('Entrada perturbación sistema')
xlabel('tiempo (minutos)')
ylabel('Perturbación (Adim)')


%-----------------------------------------------------
% Realimentación por error en estado estacionario

dif = 1e-05;
%pr = [0.85 -0.95 0.9];
pr = [0.96 0.82 0.82];
[K,L] = poles_int(sysD,pr);
sim('RealNoLineal')


%---------------------------------------------
% GSUA 

model='resistenciaBacterias_GSUA';
ParIn={'br','bs','q'};
rMethod='range';
Ranges={0 2; 0 2; 0 1};
SampleMethod='LatinHypercube';
SensMethod='Saltelli';
y_exp='';
N=500;
parallel=1;
sim_progress='on';
nominals=[0.14 0.1 4^(-3)];
additional={"Dynamic", 0.5};
[M,Par,tOut,y_Nom] = sens_Data_Prep(model,ParIn,SampleMethod,N,Ranges,rMethod,nominals,additional)

figure(1)
clf %Limpia el contenido actual de la figura
sens_plot('ScatterParameter',Par,SensMethod,M)