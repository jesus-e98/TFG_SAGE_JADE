%% LABORATORIO 5G - SWAT TIC-244
% Universidad de Granada
%--------------------------------------------------------------------------
% Algoritmo JADE (Joint Azimuth, Elevation and Delay Estimation) para
% localización 3D en interiores adaptado a las medidas de la cámara
% anecoica
%--------------------------------------------------------------------------
% Trabajo de Fin de Grado - Curso 2019-2020 -
% Realizado por: Jesús Enrique Fernández-Aparicio Ortega
% Tutor: Juan Francisco Valenzuela Valdés 
%--------------------------------------------------------------------------
% Llamada: [aml_azm,aml_ele,aml_toa,aml_amp] = jade_alg(medida,L,iteraciones,polar,paramS,toa_max)
% 
% Parámetros:
%   -medida: nombre del fichero de medida de la cámara ('nombre_fichero')
%   -L: número de MPC a estimar. 
%   -iteraciones: número de repeticiones del proceso de búsqueda
%   (recomendado [1-3])
%   -polar: seleccionamos la posición del vector de polarización deseada.
%   Se puede dejar en blanco y se tomará el valor central sin polarización.
%   -paramS: seleccionamos la dirección del canal que queremos medir
%   correspondiente a uno de los parámetros S
%           -'S21': parámetro S21
%           -'S12': parámetro S12
%   -toa_max: valor máximo del rango de búsqueda de MPCs (ns)
%
% Salida:
%   -aml_azm: ángulo azimutal phi de llegada para cada MPC estimado (grados)
%   -aml_ele: ángulo de elevación theta de llegada para cada MPC estimado (grados)
%   -aml_toa: tiempo de llegada para cada MPC estimado (s)
%   -aml_amp: amplitud de llegada para cada MPC estimado (dB)
%--------------------------------------------------------------------------

function [aml_azm,aml_ele,aml_toa,aml_amp] = jade_alg(medida,L,iteraciones,polar,paramS,toa_max)

%-------------PARA PRUEBAS-------------------------------------------------
% clear all; clc;
% tic
% medida = 'Medidas/<fichero_medida>';
% L = ;
% iteraciones = ;
% polar = ;
% paramS = ;
% toa_max = ;
% snr = ;
%--------------------------------------------------------------------------

%Cargamos la medida
load(medida)

%Seleccionamos el valor correspondiente al parámetro S
switch(paramS)
    case 'S12' %Parámetro S21
        Param_S = 1;
    case 'S21' %Parámetro S12
        Param_S = 2;
    otherwise
        Param_S = 1;
end

%Rango de frecuencias
FrecuenciaStart=freqInicial*1e9;
FrecuenciaStop=freqFinal*1e9;
freq_n = linspace(FrecuenciaStart,FrecuenciaStop,puntosFreq);

%Separación de cada frecuencia
if length(freq_n)>1
    delta_f = abs(freq_n(1) - freq_n(2));
else
    delta_f = 0;
end
lambda = 3e8./mean(freq_n);


M = size(S,2)*size(S,1); %Número de elementos
K = length(freq_n); %Número de frecuencias

%Obtenemos el array de posición de cada uno de los elementos del array de
%antenas en 3D en el receptor
array_pos = [];
for i = 1:length(Xs)
    for j = 1:length(Zs)
        aux = [Xs(i),Zs(j)]*1e-2;
        array_pos = [array_pos; aux];
    end 
end

%Fases correspondientes a cada frecuencia y a cada posición en el array
K_tof = -1i*2*pi*(freq_n).';
K_pos = -1i*2*pi/lambda*array_pos;

%Para cada una de las posiciones del transmisor en la mesa (Ys) en grados
for y = 1:length(Ys)
    %Iteraciones por defecto
    for mc=1:100
        %Tomamos la medida del canal concreto
        X = reshape(S(:,:,y,polar,:,Param_S),M,K);

        %-- Comienzo del algoritmo JADE------------------------------------
        %Inicialización vacía
        theta_l = [];  phi_l   = [];
        tof_l   = [];  beta_l  = [];
        
        %Estimación para cada MPC
        for l=1:L 
            
            if isempty(theta_l) %Inicialización en la primera iteración
                X_res = X;
            else %Aplicación de SIC en el resto de iteraciones
                for k=1:K
                    X_res(:,k) = X(:,k) - arst(K_pos,theta_l,phi_l)*(exp(K_tof(k)*tof_l).*beta_l)';
                end
            end
            
            %Obtención de la matriz de covarianza
            Rx  = X_res*X_res'/K;
            Rx2 = sqrtm(Rx);
            
            %Inicialización para estimación de DOA
            azm_range = 0:10:180;
            ele_range = 0:10:180;
            for ke=1:length(ele_range)
                spec(ke,:) = sum(abs(arst(K_pos,azm_range,ele_range(ke))'*Rx2).^2, 2);
            end
            %Maximización
            [azm, ele] = find_max_peak(spec, azm_range, ele_range);

            theta_l = [theta_l, azm];
            phi_l   = [phi_l, ele];

            %Hallamos steering vector y la matriz U
            A_l = arst(K_pos, theta_l, phi_l);
            U = (A_l'*A_l)\A_l'*X;

            %--- Estimación de TOA-----------------------------------------
            %Precisión de la búsqueda
            detaTOF = [1, 0.5, 0.1]*1e-9;
            %Rangos de valores a buscar
            rangTOF = [0, toa_max; 5, 5; 0.5, 0.5]*1e-9;
            
            %Para cada MPC que tenemos estimado
            for l2=1:size(U,1)
                tof_l(l2) = 0;
                
                %Búsqueda con rangos más estrechos y con mayor precisión en
                %cada iteración
                for step=1:length(detaTOF)
                    %Rango de valores
                    tof_range = tof_l(l2)-rangTOF(step,1):detaTOF(step):tof_l(l2)+rangTOF(step,2);
                    %Maximización
                    cost_toa = abs(exp(K_tof*tof_range)'*U(l2,:).');
                    [mv,mi] = max(cost_toa);
                    tof_l(l2) = tof_range(mi);
                end
                %Hallamos la atenuación del MPC
                beta_l(l2) = exp(K_tof*tof_l(l2))' * U(l2,:).'/K;
            end

            %--- Estimación de ángulos--------------------
            rr = exp(-1i*2*pi*freq_n.'*tof_l).' ;
            B = X*rr'/(rr*rr');
            
            %Rango de búsqueda progresivamente menor
            detaDOA = [2, 0.2, 0.02];
            
            %Para cada MPC
            for l3=1:size(B,2)
                
                %Buscamos en un rango progresivamente menor
                for step=1:length(detaDOA)
                    azm_range = theta_l(l3)-detaDOA(step)*5:detaDOA(step):theta_l(l3)+detaDOA(step)*5;
                    ele_range = phi_l(l3)-detaDOA(step)*5:detaDOA(step):phi_l(l3)+detaDOA(step)*5;
                    %Maximización
                    [theta_l(l3), phi_l(l3)] = az_el_estimation(B(:,l3), azm_range, ele_range, K_pos);
                end                    
            end
        end      
    end
    %Preparamos los datos finales
    aml_azm(y,:)=phi_l;
    aml_ele(y,:)=theta_l;
    aml_toa(y,:)=tof_l;
    aml_amp(y,:)=10.*log10(abs(beta_l));
    
end
% Fin del algoritmo
end

%------------------PARA REPRESENTACIÓN SI SE EJECUTA COMO SCRIPT-----------
% ytx = 3;
% figure(1)
% subplot(1,3,1)
% scatter3(aml_azm(ytx,:),aml_toa(ytx,:)*1e9,aml_ele(ytx,:),[],10*log10(abs(aml_amp(ytx,:))),'filled');
% xlabel("Azimuth (º)"); ylabel("TOA (ns)"); zlabel("Elevation (º)");view(0,0)
% xlim([0 180]); zlim([0 180]); 
% grid on
% colormap(jet)
% colorbar
% set(gca,'ColorScale','log');

% subplot(1,3,2)
% scatter3(aml_azm(ytx,:),aml_toa(ytx,:)*1e9,aml_ele(ytx,:),[],10*log10(abs(aml_amp(ytx,:))),'filled');
% xlabel("Azimuth (º)"); ylabel("TOA (ns)"); zlabel("Elevation (º)");view(90,90)
% xlim([0 180]); zlim([0 180]); 
% grid on
% colormap(jet)
% colorbar
% set(gca,'ColorScale','log');

% subplot(1,3,3)
% scatter3(aml_azm(ytx,:),aml_toa(ytx,:)*1e9,aml_ele(ytx,:),[],10*log10(abs(aml_amp(ytx,:))),'filled');
% xlabel("Azimuth (º)"); ylabel("TOA (ns)"); zlabel("Elevation (º)");view(90,0)
% xlim([0 180]); zlim([0 180]); 
% grid on
% colormap(jet)
% colorbar
% set(gca,'ColorScale','log');
%--------------------------------------------------------------------------
