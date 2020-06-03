%%%%%%%
%Eduardo Cabrera Lopez
%NUA:146618
%%%%%%%

clc;
clear all;
close all;
syms x t n wt t
A = [-pi pi] ;
f = sin(2*t)+ cos(3*t)-sin(t); %Funcion a ver cual se deja
x = linspace(min(A), max(A), 1000);
fx = 0;
a = 3; %Numero Armonicos


%Aqui verifico si la funcion sera par o impar .-.
for i=1:length(A)-1
    if mod(i,2) == 1
        %Impar        
        fx = fx +((x>=A(i))&(x<=A(i+1))).*subs(f(i),x);
    else
        %Par        
        fx = fx +((x>A(i))&(x<A(i+1))).*subs(f(i),x);
    end
end
subplot(2,1,1)
plot(x,fx, 'black', 'Linewidth' , 1.3);
hold on
plot(x+max(x)-min(x), fx , 'black', 'Linewidth' , 1.3);
plot(x-max(x)+min(x), fx , 'black', 'Linewidth' , 1.3);
plot([max(x) max(x)], [fx(1) fx(end)], 'black', 'Linewidth' , 1.3);
plot([min(x) min(x)], [fx(end) fx(1)], 'black', 'Linewidth' , 1.3);
grid on;
xlabel('Tiempo');
ylabel('Amplitud');
title(['Funcion Origial: $' latex(f) '$'],'Interpreter','latex')
pause(0.00001)

f = sym(f);
T = max(A)-min(A);
wo = 2*pi/(T);

Ao = 0;
for i = 1:length(f)
     Ao = Ao + int(f(i),'t',A(i),A(i+1));
     disp('Ya termine Ao') %METODO DEBUG ANTI CIERRES 
end

An = 0;
for i = 1:length(f)
    An = An + int(f(i)*cos(n*wo*t),'t',A(i),A(i+1));
    disp('Ya termine An') %METODO DEBUG ANTI CIERRES 
end

Bn = 0;
for i = 1:length(f)
    Bn = Bn + int(f(i)*sin(n*wo*t),'t',A(i),A(i+1));
    disp('Ya termine Bn') %METODO DEBUG ANTI CIERRES
    
end

t = linspace(min(A)-T, max(A)+T, 1000);
ft = zeros(a,1000);

for i=1:a    
    ft(i,:) = (subs(Bn, 'n', i).*sin(i*wo*t))+(subs(An,'n',i).*cos(i*wo*t));
    subplot(2,1,2)
    plot(t , ft(i,:), 'Linewidth' , 1.3)

    hold on
    grid on
end
title("Total de armonicos "+ num2str(a));
disp('Aun no me congelo... Aun sigo calculando') %METODO DEBUG ANTI CIERRES 
saveas(gcf,'Funcion.pdf')
ft = zeros(a,1000);
figure
for i=1:a  
    ft(i,:) = (subs(Bn, 'n', i).*sin(i*wo*t))+(subs(An,'n',i).*cos(i*wo*t));
    subplot(a,1,i)
    plot(t, Ao+sum(ft), 'Linewidth' , 1.3)
    title("Suma de "+num2str(i)+" Armonico(s)")    
 
    hold on
    grid on;
end
saveas(gcf,'Armonicos.pdf')


% plot(t, Ao+sum(ft), 'Linewidth' , 1.3)
% grid on;
% xlabel('Tiempo');
% ylabel('Amplitud');
% title('Sumatoria de Armonicos de Fourier')



    
