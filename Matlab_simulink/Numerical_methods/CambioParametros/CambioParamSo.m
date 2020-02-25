clear all
% Este codigo saca la grafica, y la tabla de los datos y la escribe en 
% una tabla de excel llamada filename
filename = 'CambioParamSoSalidas.xlsx';
tiempoFinal = 80;
numSim = 4;
n = 1602;
load_system('resistenciaBacterias');
Y = zeros(n,numSim);
So = 0:0.5:1.5;
str = '';
for i = 0:(numSim-1)
    set_param('resistenciaBacterias/Subsystem','par','[0.14,0.1,0.3960,4^(-3),0.2,0.09,0.0083]','x0',strcat('[',num2str(So(i+1)),',0,0]'));
    set_param('resistenciaBacterias','StopTime',num2str(tiempoFinal));
    sim('resistenciaBacterias');
    str = strcat(str,',So = ',num2str(So(i+1)));
    Y1(:,i+1) = Rout;
    Y2(:,i+1) = Rout + Sout;
end
% plot(tout,Y1)
% title ('Comparación bacterias resistentes cambio de So')
% ylabel('Bacterias resistentes')
plot(tout,Y2)
title ('Comparación bacterias sensibles + bacterias resistentes cambio de So')
ylabel('Bacterias bacterias sensibles + bacterias resistentes')

xlabel('Tiempo en minutos')
legend('So = 0','So = 0.5','So = 1','So = 1.5') 
A = [So',mean(Y1)',mean(Y1(n-100:n,:))',min(Y1)',max(Y1)',mean(Y2)',mean(Y2(n-100:n,:))',min(Y2)',max(Y2)'];
A = round(A,6)
T = array2table(A,...
    'VariableNames',{'ValorSo','ValorMediaY1','Valor_Estab_Y1',...
    'minValorY1','maxValorY1','ValorMediaY2','Valor_Estab_Y2',...
    'minValorY2','maxValorY2'})
writetable(T,filename,'Sheet',1,'Range','A1')