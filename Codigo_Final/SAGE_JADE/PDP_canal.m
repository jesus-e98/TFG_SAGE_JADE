%% LABORATORIO 5G - SWAT TIC-244
% Universidad de Granada
%--------------------------------------------------------------------------
% Función para hallar el PDP medio recibido en una dirección y polarización
% determinada. Obtenemos un par de gráficas donde se muestra el PDP en
% unidades lineales y logarítmicas.
%--------------------------------------------------------------------------
% Trabajo de Fin de Grado - Curso 2019-2020 -
% Realizado por: Jesús Enrique Fernández-Aparicio Ortega
% Tutor: Juan Francisco Valenzuela Valdés 
%--------------------------------------------------------------------------
% Llamada: [RMS_DS, RMS_DS_Envolvente] = PDP_canal(medida,pos_y,polar,canal)
%
% Parámetros:
%   -medida: nombre del fichero de medida de la cámara ('nombre_fichero')
%   -pos_y: índice de una posición del Eje Y de la cámara. 
%   -polar: seleccionamos la posición del vector de polarización eje E deseada.
%   Se puede dejar en blanco y se tomará el valor central sin polarización.
%   -canal: seleccionamos la dirección del canal que queremos medir
%   correspondiente a uno de los parámetros S
%           -'S21': parámetro S21
%           -'S12': parámetro S12
%
% Salida:
%   -RMS_DS: valor medio de retardo del PDP
%   -RMS_DS_Envolvente: valor medio de retardo de la envolvente de PDP
%--------------------------------------------------------------------------

function [RMS_DS, RMS_DS_Envolvente] = PDP_canal(medida,pos_y,polar,canal)

%% CARGA DE DATOS CÁMARA %%

switch(canal)
    case 'S12'
        param = 1;
    case 'S21'
        param = 2;
    otherwise
        param = 1;
end

load(medida)
Sxy=S;
total = zeros(1,size(Sxy,5));


%Hallamos el PDP medio de todos los elementos de la antena
for i = 1:length(Xs)
    for j = 1:length(Zs)
        p_s12 = reshape(Sxy(i,j,pos_y,polar,:,param),1,size(Sxy,5));
%         p_s12 = awgn(p_s12,40,'measured');
        total(1,:) = total(1,:) + p_s12(1,:); 
    end
end
p_s12 = total./(length(Xs)*length(Zs))

%Rango de frecuencias
fini=freqInicial*1e9;
fstop=freqFinal*1e9;
fres_aux=linspace(fini,fstop,puntosFreq);


%% Relleno con ceros y reflejo de la señal para que sea simétrica %%
dF=fres_aux(2)-fres_aux(1); %Resolución en frecuencia
T=1/(2*fstop);  % Periodo de Muestreo
N_frec=round(fstop/dF);  % Desde 0
relleno=zeros(1,N_frec-length(p_s12));
P_rellena=[relleno p_s12];

Preflejada(1:length(P_rellena))=P_rellena(end:-1:1);



%% IFFT %%
[h_t3] = ifft([P_rellena conj(Preflejada)],'symmetric');
tiempo=1:length(h_t3);
t3=T*tiempo;

h_t=h_t3; t=t3;

%% Calculo del PDP y Normalización %%
PDP=abs(h_t).^2;
PDP=PDP(1:round(length(PDP)/2));    %Por simetría de la IFFT, me quedo con la primera mitad
[max_p ind_max] =max(abs(PDP));     %Normalización
%PDP=PDP(ind_max:end);
PDP_norm=abs(PDP)/max_p;
[u,b]=max(PDP_norm);
t=t(1:length(PDP_norm));


%% Representación
figure(10)
subplot(2,1,1)
plot(t/1e-9,PDP_norm);hold on
title('Power Delay Profile (PDP)');
xlabel('Tiempo (ns)');
ylabel('Potencia Normalizada');
%xlim([0 500]); 
%ylim([0 1]); 

%% Señal en dB %%

PDP_log=10*log10(PDP_norm);
figure(10)
subplot(2,1,2)
hold on
plot(t/1e-9,PDP_log);hold on
title('Power Delay Profile (PDP)');
xlabel('Tiempo (ns)');
ylabel('Potencia Normalizada (dB)');
ylim([-80 0]); 

%% CÁLCULO DEL RMS DELAY SPREAD A PARTIR DE LA SEÑAL OBTENIDA CON IFFT %%
%ds_bluetest=sqrt((sum(PDP_norm.*(t.^2))/sum(PDP_norm))-(sum(PDP_norm.*t)/sum(PDP_norm)).^2)
tau_med=sum(PDP_norm.*t)/sum(PDP_norm);
RMS_DS=sqrt(sum(PDP_norm.*((t-tau_med).^2))/sum(PDP_norm))

%% Suavizado de la señal (Me quedo con la envolvente) %%
PDP_log = env_secant(t, 10*log10(PDP_norm),round(length(P_rellena)*0.01),'top');
plot(t/1e-9,PDP_log,'r');hold off

%% Reconstruccion del PDP a partir de la señal suavizada %%
figure(10)
subplot(2,1,1)
PDP_log=PDP_log(1:length(t));
h_lineal=10.^(PDP_log./10);
plot(t/1e-9,h_lineal,'r');hold off
PDP_norm=h_lineal;

%% CÁLCULO DEL RMS DELAY SPREAD A PARTIR DE LA ENVOLVENTE %%
%ds_bluetest2=sqrt((sum(PDP_norm.*(t.^2))/sum(PDP_norm))-(sum(PDP_norm.*t)/sum(PDP_norm)).^2)
tau_med=sum(PDP_norm.*t)/sum(PDP_norm);
RMS_DS_Envolvente=sqrt(sum(PDP_norm.*((t-tau_med).^2))/sum(PDP_norm))

end