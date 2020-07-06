%% LABORATORIO 5G - SWAT TIC-244
% Universidad de Granada
%--------------------------------------------------------------------------
% Función para transformar de coordenadas cartesianas a esféricas y
% viveversa
%--------------------------------------------------------------------------
% Trabajo de Fin de Grado - Curso 2019-2020 -
% Realizado por: Jesús Enrique Fernández-Aparicio Ortega
% Tutor: Juan Francisco Valenzuela Valdés 
%--------------------------------------------------------------------------
% Llamada: [az_out,el_out,s_out] = az_el_solid_angles(<az_in>,<el_in>,<s_in>,modo)
% 
% Parámetros:
%   -az_in: ángulo azimutal phi (grados)
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


function [az_out,el_out,s_out] = az_el_solid_angles(az_in,el_in,s_in,modo)
 
%Pasamos a radianes los ángulos
az_in = az_in .* pi/180;
el_in = el_in .* pi/180;

switch(modo)
    case 'r' %Cambiamos de ángulo sólido a az y el
        el_out = acos(s_in(3))*(180/pi);
        az_out = asin(s_in(2)/sin(el_out))*(180/pi);
        s_out = s_in;
    case 'd' %Cambiamos de az y el a ángulo sólido
        s_out(:) = [cos(az_in).*sin(el_in),sin(az_in).*sin(el_in),cos(el_in)]';
        az_out = az_in.*(180/pi);
        el_out = el_in.*(180/pi);
end
        

    
