%% LABORATORIO 5G - SWAT TIC-244
% Universidad de Granada
%--------------------------------------------------------------------------
% Función auxiliar del algoritmo JADE para hallar los ángulos de llegada
% mediante la maximización de la matriz B
%--------------------------------------------------------------------------
% Trabajo de Fin de Grado - Curso 2019-2020 -
% Realizado por: Jesús Enrique Fernández-Aparicio Ortega
% Tutor: Juan Francisco Valenzuela Valdés 
%--------------------------------------------------------------------------
% Llamada: [azm, ele] = az_el_estimation(b, azm_range, ele_range, array_pos)
% 
% Parámetros:
%   -b: matriz B 
%   -el_in: ángulo de elevación theta (grados)
%   -s_in: array con las posiciones en sistema cartesiano [X,Y,Z]
%   -modo: modo de operacion de la función
%           -'r': cambio de cartesianas a ángulos de esféricas
%           -'d': cambio de angulos de esféricas a cartesianas
% Salida:
%   -az_out: ángulo azimutal phi de salida correspondiente a unas coordenadas
%   cartesianas (grados)
%   -el_out: ángulo de elevación theta de salida correspondiente a unas 
%   coordenadas cartesianas (grados)
%   -s_out: array con las posiciones en sistema cartesiano correspondientes
%   a unos ángulos de esféricas [X,Y,Z]
%--------------------------------------------------------------------------

function [azm, ele] = az_el_estimation(b, azm_range, ele_range, array_pos)
    
d2pi = pi/180;
    
for ele=1:length(ele_range)
    for azm = 1:length(azm_range)
        %-----NOTA---------------------------------------------------------
        % Otras configuraciones de rho son posibles en distintos planos. La
        % expresión general sería
        % rho = [cos(azm_range*d2pi).*cos(ele_range(ke)*d2pi); sin(azm_range*d2pi).*cos(ele_range(ke)*d2pi); cos(azm_range(azm)*0).* cos(ele_range(ele)*d2pi)];
        %------------------------------------------------------------------
        rho = [cos(azm_range(azm)*d2pi).*sin(ele_range(ele)*d2pi); cos(azm_range(azm)*0).* cos(ele_range(ele)*d2pi)];
        azel(azm,ele) = pow2(abs(b'*exp(array_pos*rho)));
    end
end

%Buscamos en 2D el valor máximo
[m1, im1] = max(azel);
[m2, im2] = max(m1);
azm = azm_range(im2);
ele = ele_range(im1(im2));

%-----------------EXTRA----------------------------------------------------
% Representación en 3D de la función que estamos maximizando
% [AZ,EL] = meshgrid(azm_range, ele_range);
% figure(10); mesh(AZ,EL,azel/max(azel(:))); 
% xlabel('Azimuth (º)'); 
% ylabel('Elevation (º)');
% ylim([azm_range(1), azm_range(end)]); zlim([0 1.1]);%view(90,0);
%--------------------------------------------------------------------------
end