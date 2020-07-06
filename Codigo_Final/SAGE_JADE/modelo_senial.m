%% LABORATORIO 5G - SWAT TIC-244
% Universidad de Granada
%--------------------------------------------------------------------------
% Implementa un modelo de señal usado en el contexto del proceso
% de optimización de SAGE. Se puede cambiar por otro modelo si es
% necesario o introducir nuevos
%--------------------------------------------------------------------------
% Trabajo de Fin de Grado - Curso 2019-2020 -
% Realizado por: Jesús Enrique Fernández-Aparicio Ortega
% Tutor: Juan Francisco Valenzuela Valdés 
%--------------------------------------------------------------------------
% Llamada: x_out = modelo_senial(Theta,Dim,freq_n,array_pos_rx,modelo)
% 
% Parámetros:
%   -Theta: matriz de parámetros
%   -Dim: dimensión de la señal
%   -freq_n: vector de frecuencias
%   -array_pos_rx: posiciones de los elementos del array receptor en
%   cartesianas
%   -modelo: modelo de señal a utilizar 
%
% Salida:
%   -x_out: señal sintetizada a partir de los parámetros
%--------------------------------------------------------------------------


function x_out = modelo_senial(Theta,Dim,freq_n,array_pos_rx,modelo)


switch(modelo)
    case 'modelo1'
        %Extraemos los parámetros
        tau = Theta(1);
        azimuth = Theta(2);
        ampl = Theta(3);

        %Hallamos el steering vector
        arst_comp = [1:Dim];

        %Suponemos 0.5 lambda de separacion
        for indice_Dim = 1:length(arst_comp)
            A = exp(-1i * 2 * pi * (arst_comp(indice_Dim)-1) * 0.5 * sin(azimuth));
        end

        %Aplicamos el modelo de la señal
        x_out = A * ampl * exp(-1i * 2 * pi .* freq_n * tau);
        
    case 'modelo2'
        array_pos_tx = [0 0 0]; %Esto estaría en los parámetros en un futuro
        
        %Extraemos los parámetros
         tau = Theta(1);
         az_rx = Theta(2);
         el_rx = Theta(3);
         az_tx = Theta(4);
         el_tx = Theta(5);
         ampl = Theta (6);
         
         [~,~,eRX] = az_el_solid_angles(az_rx,el_rx,[],'d');
         f_rad_RX = diagrama_rad(az_rx,el_rx,'bocina');
         [~,~,eTX] = az_el_solid_angles(az_tx,el_tx,[],'d');
         f_rad_TX = diagrama_rad(az_tx,el_tx,'bocina');
         
         x_out = zeros(length(array_pos_tx(:,1)),length(array_pos_rx(:,1)),length(freq_n),1);
         
         for fq = 1:length(freq_n)
             %Calculamos el steering vector del transmisor
             c1 = [];
             for i = 1:length(array_pos_tx(:,1))
%                 c1_aux = f_rad_TX*exp(1i*2*pi*(freq_n(fq)/(3e8))*dot(array_pos_tx(i,:),eTX));
                c1_aux =f_rad_TX* exp(1i*2*pi*(freq_n(fq)/(3e8))*dot(array_pos_tx(i,:),eTX));
                c1 = [c1; c1_aux];
             end
             
             %Calculamos el steering vector del receptor
             c2 = [];
             for i = 1:length(array_pos_rx(:,1))
%                 c2_aux = f_rad_RX*exp(1i*2*pi*(freq_n(fq)/(3e8))*dot(array_pos_rx(i,:),eRX));
                c2_aux =f_rad_RX* exp(1i*2*pi*(freq_n(fq)/(3e8))*dot(array_pos_rx(i,:),eRX));
                c2 = [c2; c2_aux];
             end
             
             x_out(:,:,fq) = ampl*(c1*c2')*exp(-1i*2*pi*freq_n(fq)*tau);
         end

end