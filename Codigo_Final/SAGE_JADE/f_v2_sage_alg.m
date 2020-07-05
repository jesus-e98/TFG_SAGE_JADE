%% LABORATORIO 5G - SWAT TIC-244
% Universidad de Granada
%--------------------------------------------------------------------------
% Algoritmo SAGE (Space Alternating Generalized Expectation-maximization) para
% localización 3D 
%--------------------------------------------------------------------------
% Trabajo de Fin de Grado - Curso 2019-2020 -
% Realizado por: Jesús Enrique Fernández-Aparicio Ortega
% Tutor: Juan Francisco Valenzuela Valdés 
%--------------------------------------------------------------------------
% Llamada: [sage_aoa,sage_eoa,sage_toa,sage_amp] = f_v2_sage_alg(Medidas,L,iters,polar,Canal,toa_max)
% 
% Parámetros:
%   -Medidas: nombre del fichero de medida de la cámara ('nombre_fichero')
%   -L: número de MPC a estimar. 
%   -iters: número de repeticiones del proceso de búsqueda
%   (recomendado [1-3])
%   -polar: seleccionamos la posición del vector de polarización deseada.
%   Se puede dejar en blanco y se tomará el valor central sin polarización.
%   -Canal: seleccionamos la dirección del canal que queremos medir
%   correspondiente a uno de los parámetros S
%           -'S21': parámetro S21
%           -'S12': parámetro S12
%   -toa_max: valor máximo del rango de búsqueda de MPCs (ns)
%
% Salida:
%   -sage_aoa: ángulo azimutal phi de llegada para cada MPC estimado (grados)
%   -sage_eoa: ángulo de elevación theta de llegada para cada MPC estimado (grados)
%   -sage_toa: tiempo de llegada para cada MPC estimado (s)
%   -sage_amp: amplitud de llegada para cada MPC estimado (dB)
%--------------------------------------------------------------------------

function [sage_aoa,sage_eoa,sage_toa,sage_amp] = f_v2_sage_alg(Medidas,L,iters,polar,Canal,toa_max)

%--------PARA REALIZACIÓN DE PRUEBAS --------------------------------------
% Comentar la línea function y end.
% Descomentar las líneas a continuación
% clc; clear all;
% Medidas = 'Medidas/<nombrefichero>';
% L = ;
% iters = ;
% polar = ;
% SNR = ;
% Canal = ;
% toa_max = ;
%--------------------------------------------------------------------------

load(Medidas);
switch(Canal)
    case 'S12' %Parámetro S12
        Param_S = 1;
    case 'S21' %Parámetro S21
        Param_S = 2;
    otherwise
        Param_S = 1;
end
freq_n = linspace(freqInicial,freqFinal,puntosFreq)*1e9;
lambda = 3e8/mean(freq_n);

%Dimensiones de la matriz de medida
K = length(freq_n);
M = length(Xs)*length(Zs);   

%Obtenemos el array de posición de cada uno de los elementos del array de
%antenas en 3D en el receptor
array_pos_rx = [];
for i = 1:length(Xs)
    for j = 1:length(Zs)
        aux = [Xs(i),0,Zs(j)]*1e-2;
        array_pos_rx = [array_pos_rx; aux];
    end 
end

%Hacemos lo mismo para un array de transmisión
%*************************FUTURA AMPLIACIÓN********************************
%Más adelante podremos hacer lo mismo si tenemos un RPA en el emisor por
%ahora se queda así
array_pos_tx = [0 0 0];
%**************************************************************************


%Para el estudio de las medidas se usará un doble bucle que obtendrá los
%MPCs para valores de  
% -Orientación/ángulo del Transmisor (Ys)
% -Ruido que añadiremos por software (SNR)
%**************POSIBLE MEJORA**********************************************
% -Distintos Valores de Polarización (Es) -> Está fijada a 0º de desvío
%**************************************************************************

for y=1:length(Ys)
    %Tomamos la señal a utilizar
    signal_X = reshape(S(:,:,y,polar,:,Param_S),M,K); 
%     signal_X = awgn(signal_X,SNR);
    
    %Inicializamos a cero los valores a estimar
    toa = zeros(1,L); 
    OmegaR = zeros(2,L); %Estos parámetros contienen a su vez los ángulos phi y
    OmegaT = zeros(2,L); %theta azimutal y de elevacion respectivamente
    ampl = zeros(1,L);
    
    for iteracion = 1:iters
        fprintf("\nIteracion: %i/%i\n",iteracion,iters);
        if iteracion == 1   %Inicialización del algoritmo
            for IC = 1:L
                %Cancelación de componentes usando técnica SIC
                x_i = intercancelacion(signal_X,toa,OmegaR,OmegaT,ampl,IC,freq_n,array_pos_rx,'serie');
                
                %Resolución de búsqueda -> cte
                deltaTOF = [0.5, 0.5, 0.5]*1e-9;
                
                %Rango de valores a buscar
                rangTOF=[   0, toa_max;      %step1
                            5, 5;       %step2
                            1, 1]*1e-9; %step3
                
                for step = 1:length(deltaTOF)
                    toa_range = toa(IC)-rangTOF(step,1):deltaTOF(step):toa(IC)+rangTOF(step,2);
                    %Filtramos los valores que no sean positivos
                    toa_range = toa_range(toa_range>0);
                    
                    %Cálculo de la función de coste
                    Coste = zeros(1,length(toa_range));
                    for t = 1:length(Coste)
                        Coste(t) = sum(pow2(abs(sum(exp(1i*2*pi*freq_n*toa_range(t)).*permute(x_i,[2 1])',2))),1);
                    end
                    %Maximización
                    [~,Index] = max(Coste);
                    toa(IC) = toa_range(Index);
                end
                
                %Búsqueda de DOA con mayor precisión
                deltaDOA = [10, 2, 0.2, 0.02];
                
                for step = 1:length(deltaDOA)
                    %Definimos los rangos a buscar
                    
                    if deltaDOA(step) == 10 %Si estamos en la primera iteración
                        azm_range = 0:10:180;
                        ele_range = 0:10:180;
                    else %En el resto de iteraciones
                        azm_range = OmegaR(1,IC)-deltaDOA(step)*5:deltaDOA(step):OmegaR(1,IC)+deltaDOA(step)*5;
                        ele_range = OmegaR(2,IC)-deltaDOA(step)*5:deltaDOA(step):OmegaR(2,IC)+deltaDOA(step)*5;
                    end
                    
                    %Cálculo de la función de coste
                    Coste2 = zeros(length(azm_range),length(ele_range));
                
                    for azm = 1:size(Coste2,1)
                        for ele = 1:size(Coste2,2)
                            c2 = steering(array_pos_rx,azm_range(azm),ele_range(ele),lambda,'tx_rx');
                            Coste2(azm,ele) = pow2(abs(c2' * sum(exp(1i*2*pi*freq_n*toa(IC)).*permute(x_i,[2 1])',2)));
                        end
                    end
                    %Maximización
                    [azm_max,ele_max] = find_max_peak(Coste2,azm_range,ele_range);
                    OmegaR(1,IC) = azm_max;
                    OmegaR(2,IC) = ele_max;
                end
            end    
            %FIN DEL CICLO DE INICIALIZACIÓN
            
        else                %Parte iterativa del algoritmo
            for IC = 1:L
                 x_i = intercancelacion(signal_X,toa,OmegaR,OmegaT,ampl,IC,freq_n,array_pos_rx,'serie');
                deltaTOF = [0.5, 0.5, 0.5]*1e-9;
                
                rangTOF=[   0, toa_max;      %step1
                            5, 5;       %step2
                            1, 1]*1e-9; %step3
                
                for step = 1:length(deltaTOF)
                    if step == 1
                        toa_range = [0:0.5:toa_max]*1e-9;
                    else 
                        toa_range = toa(IC)-rangTOF(step,1):deltaTOF(step):toa(IC)+rangTOF(step,2);
                    end
                    %Filtramos los valores que no sean positivos
                    toa_range = toa_range(toa_range>0);
                    
                    Coste = zeros(1,length(toa_range));
                    for t = 1:length(Coste)
                        Theta_l = [toa_range(t) OmegaT(1,IC) OmegaT(2,IC) OmegaR(1,IC) OmegaR(2,IC)];
                        Coste(t) = funcion_coste2(x_i,freq_n,Theta_l,array_pos_rx,array_pos_tx,'3D_AOA2');
                    end
                    [~,Index] = max(abs(Coste));
                    toa(IC) = toa_range(Index);
                end
                
                deltaDOA = [10, 2, 0.2, 0.02];
                
                for step = 1:length(deltaDOA)
                    %Definimos los rangos a buscar
                    
                    if deltaDOA(step) == 10 %Si estamos en la primera iteración
                        azm_range = 0:10:180;
                        ele_range = 0:10:180;
                    else %En el resto de iteraciones
                        azm_range = OmegaR(1,IC)-deltaDOA(step)*5:deltaDOA(step):OmegaR(1,IC)+deltaDOA(step)*5;
                        ele_range = OmegaR(2,IC)-deltaDOA(step)*5:deltaDOA(step):OmegaR(2,IC)+deltaDOA(step)*5;
                    end
       
                    Coste2 = zeros(length(azm_range),length(ele_range));
                
                    for azm = 1:size(Coste2,1)
                        for ele = 1:size(Coste2,2)
                            Theta_l = [toa(IC) OmegaT(1,IC) OmegaT(2,IC) azm_range(azm) ele_range(ele)];
                            Coste2(azm,ele) = funcion_coste2(x_i,freq_n,Theta_l,array_pos_rx,array_pos_tx,'3D_AOA2');
                        end
                    end
                    [azm_max,ele_max] = find_max_peak(Coste2,azm_range,ele_range);
                    OmegaR(1,IC) = azm_max;
                    OmegaR(2,IC) = ele_max;
                end
               
                
                Theta_l = [toa(IC) OmegaT(1,IC) OmegaT(2,IC) OmegaR(1,IC) OmegaR(2,IC)];
                Coste3 = funcion_coste2(x_i,freq_n,Theta_l,array_pos_rx,array_pos_tx,'3D_AOA2');
                ampl(IC) = (1/(M*K))*Coste3;    
            end
        end
        %Almacenamos los datos de salida
        sage_amp(y,:)=ampl;
        sage_toa(y,:)=toa;
        sage_aoa(y,:)=OmegaR(1,:);
        sage_eoa(y,:)=OmegaR(2,:);
    end
    
    
    
end
toc
%Fin del algoritmo
end
%--------------PARA REPRESENTACIÓN GRÁFICA AL EJECUTAR COMO SCRIPT---------
% ytx = 1;
% figure(1)
% subplot(1,3,1)
% scatter3(sage_aoa(ytx,:),sage_toa(ytx,:)*1e9,sage_eoa(ytx,:),[],10*log10(abs(sage_amp(ytx,:))),'filled');
% xlabel("Azimuth (º)"); ylabel("TOA (ns)"); zlabel("Elevation (º)");view(0,0)
% xlim([80 100]); zlim([80 100]); 
% grid on
% colormap(jet)
% colorbar
% set(gca,'ColorScale','log');
% 
% subplot(1,3,2)
% scatter3(sage_aoa(ytx,:),sage_toa(ytx,:)*1e9,sage_eoa(ytx,:),[],10*log10(abs(sage_amp(ytx,:))),'filled');
% xlabel("Azimuth (º)"); ylabel("TOA (ns)"); zlabel("Elevation (º)");view(90,90)
% xlim([80 100]); zlim([80 100]); 
% grid on
% colormap(jet)
% colorbar
% set(gca,'ColorScale','log');
% 
% subplot(1,3,3)
% scatter3(sage_aoa(ytx,:),sage_toa(ytx,:)*1e9,sage_eoa(ytx,:),[],10*log10(abs(sage_amp(ytx,:))),'filled');
% xlabel("Azimuth (º)"); ylabel("TOA (ns)"); zlabel("Elevation (º)");view(90,0)
% xlim([80 100]); zlim([80 100]); 
% grid on
% colormap(jet)
% colorbar
% set(gca,'ColorScale','log');

