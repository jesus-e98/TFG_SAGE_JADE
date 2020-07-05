%Fundamentos
%ULA
%H(wt)
f=5e10; %Ajusta la frecuencia para verlo mejor
M=[8,64,128]; %Puedes cambiar el número de elementos
w=2*pi*f; 
Theta=[-4*pi:pi/16:4*pi]; 
%Estos datos no se usan. Asumimos separación de 0.5lambda
c=3e8;
d=c/f*0.5;
%--------------------------------------------------------
H=zeros(length(Theta),1);
figure;
title('Función de respuesta del array para diferentes valores de M');
for n=1:length(M)
    for x=1:length(Theta)
        for j=0:M(n)-1
            H(x)=H(x)+exp(1i*w*(j*0.5).*sin(Theta(x)));
        end
    end
    H=abs(fft(H/M(n)));
    plot(Theta,H);
    hold on
end
hold off
title('Función de respuesta de un ULA para diferentes valores de M');
xlabel('DOA')
ylabel('Magnitud')
legend('M=8','M=64','M=128')
xlim([-10 10])
grid on

