%% LABORATORIO 5G - SWAT TIC-244
% Universidad de Granada
%--------------------------------------------------------------------------
% Función auxiliar para el cálculo del máximo de una matriz de dos
% dimensiones.
%--------------------------------------------------------------------------
% Trabajo de Fin de Grado - Curso 2019-2020 -
% Realizado por: Jesús Enrique Fernández-Aparicio Ortega
% Tutor: Juan Francisco Valenzuela Valdés 
%--------------------------------------------------------------------------
% Llamada: [phi, theta] = find_max_peak(spec, azm_range, ele_range)
% 
% Parámetros:
%   -spec: función 2D con los valores en el eje Z para cada posición X,Y
%   -azm_range: rango de ángulos azimutales a buscar
%   -ele_range: rango de ángulos de elevación a buscar
%
% Salida:
%   -phi: ángulo azimutal phi máximo (grados)
%   -theta: ángulo de elevación theta máximo (grados)
%--------------------------------------------------------------------------

function [phi, theta] = find_max_peak(spec, azm_range, ele_range)
    
[m1, im1] = max(spec);
[m2, im2] = max(m1);
phi = azm_range(im2);
theta = ele_range(im1(im2));

%-----------PARA REPRESENTACIÓN DE LA FUNCIÓN 2D A MAXIMIZAR---------------
% [AZ,EL] = meshgrid(azm_range, ele_range);
% figure(10); mesh(AZ,EL,spec); xlabel('Azimuth (degree)'); ylabel('Elevation (degree)');
% ylim([azm_range(1), azm_range(end)]); zlim([0 1.1]);%view(90,0);

end