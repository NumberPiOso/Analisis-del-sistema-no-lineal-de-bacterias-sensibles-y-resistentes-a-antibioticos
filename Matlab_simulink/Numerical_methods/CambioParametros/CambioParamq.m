clear all
% Este codigo saca la grafica, y la tabla de los datos y la escribe en 
% una tabla de excel llamada filename
filename = 'CambioParamQ.xlsx';
tiempoFinal = 1000;
numSim = 5;
n = 20002;
load_system('resistenciaBacterias');
Y = zeros(n,numSim);
vq = 0:4^(-4):4^(-3);
str = '';
for i = 0:(numSim-1)
    set_param('resistenciaBacterias/Subsystem','par',strcat('[0.14,0.1,0.3960,',num2str(i*4^(-4)),',0.2,0.09,0.0083]'),'x0','[0.2,0 0]');
    set_param('resistenciaBacterias','StopTime',num2str(tiempoFinal));
    sim('resistenciaBacterias');
    str = strcat(str,',q = ',num2str(i*4^(-4)));
    Y(:,i+1) = Rout;
end
plot(tout,Y)
title ('Comparación bacterias resistentes cambio del parámetro q')
ylabel('Bacterias resistentes')
xlabel('Tiempo en minutos')
legend('q = 0','q = 0.0039063','q = 0.0078125','q = 0.011719','q = 0.015625') 
A = [vq',mean(Y)',mean(Y(n-100:n,:))',min(Y)',max(Y)'];
A = round(A,5)
T = array2table(A,...
    'VariableNames',{'Valorq','Valor_Promedio_R','Valor_Estab_R','minValorR','maxValorR'})
writetable(T,filename,'Sheet',1,'Range','A1')