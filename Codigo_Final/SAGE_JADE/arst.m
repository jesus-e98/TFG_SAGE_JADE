%% LABORATORIO 5G - SWAT TIC-244
% Universidad de Granada
%--------------------------------------------------------------------------
% Función para estimar la función de respuesta de un array en el plano XZ 
% usando la fórmula simple en función de una matriz de posiciones y los
% ángulos de llegada.
%--------------------------------------------------------------------------
% Trabajo de Fin de Grado - Curso 2019-2020 -
% Realizado por: Jesús Enrique Fernández-Aparicio Ortega
% Tutor: Juan Francisco Valenzuela Valdés 
%--------------------------------------------------------------------------
% Llamada: a = arst(K_pos, azm, ele)
% 
% Parámetros:
%   -K_pos: matriz con las posiciones en coordenadas cartesianas en 
%    2D de los elementos del array.
%   -azm: array de ángulos azimutales a estudiar su respuesta
%   -ele: array de ángulos de elevación a estudiar su respuesta
% 
% Salida:
%   -a: respuesta del array de antenas a los ángulos azm y ele en cada
%   elemento
%--------------------------------------------------------------------------

function a = arst(K_pos, azm, ele)

d2pi = pi/180; %Transformación de ángulo en grados a radianes

%Consideramos el plano XZ
rho=[cos(azm.*d2pi).*sin(ele.*d2pi); cos(0.*azm).*cos(ele.*d2pi)];

%--------------EXTRA-------------------------------------------------------
% Considerando coordenadas en 3D tendríamos
% rho = [cos(azm*d2pi).*cos(ele*d2pi); sin(azm*d2pi).*cos(ele*d2pi); cos(0.*azm).*cos(ele.*d2pi)];
%--------------------------------------------------------------------------

a = exp(K_pos*rho); %Obtenemos la respuesta del array

end