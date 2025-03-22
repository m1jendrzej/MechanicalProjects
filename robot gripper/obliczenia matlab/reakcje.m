clc; clear;
L=40;
z_plot=30:0.01:52.5; %%%%Zakres wymiarów chwytanego przedmiotu dzielony na 2
N=[0 175]';
a=15;
b=22;
d=79;
c=10;
m=8.0;

omega =[0 -1; 
        1 0];

vector_c=[0 10]';
vector_r3=[0 -(b-c)]';

z_plot=2*z_plot;
i=1;
for z=30.0:0.01:52.5
    
    rn=52.2+m;

    alfa = asin((d-(z+c))/L);
    beta = -asin(((z+c)-(a+b))/L);
    
    vector_rn=[rn 0]';
    
     b_vector=[-N; -[omega*vector_rn]'*N];
     
     A=[cos(alfa) cos(beta) cos(beta);
        sin(alfa) sin(beta) sin(beta);
        (omega*vector_c)'*[cos(alfa) sin(alfa)]' (omega*vector_r3)'*[cos(beta) sin(beta)]' (omega*vector_rn)'*[cos(beta) sin(beta)]'];
     
     
     q=A\b_vector;
%###########################################################################################
%DANE DO WYKRESÓW
%###########################################################################################

%siły reakcji

     R_1(i)=q(1);
     R_2(i)=q(2);
     R_3(i)=q(3);
     
     R_1x(i)=R_1(i)*cos(alfa);
     R_1y(i)=R_1(i)*sin(alfa);

     R_2x(i)=R_2(i)*cos(beta);
     R_2y(i)=R_2(i)*sin(beta);
     
     R_3x(i)=R_3(i)*cos(beta);
     R_3y(i)=R_3(i)*sin(beta);
    
    %siła na siłownik
    S(i)=2*R_1(i)*cos(alfa);


    i=i+1;
end

i
R_1

%###########################################################################################
%Reakcja R_1
%###########################################################################################
% 2. Tworzenie wykresu
figure; % Otwiera nowe okno wykresu
plot(z_plot, R_1, '-g', 'LineWidth', 3); 
hold on; 
plot(z_plot, R_1x, '-r', 'LineWidth', 1); 
hold on;
plot(z_plot, R_1y, '-b', 'LineWidth', 1);
hold off;

% 3. Dodanie opisu osi i tytułu
xlabel('szerokosc przedmiotu z [mm]'); % Opis osi X
ylabel('Siła Reakcji R1 i jej składowe [N]'); % Opis osi Y
title('Siła reakcji R_1'); % Tytuł wykresu

% 4. Dodanie legendy
legend('R_1', 'R_1x', 'R_1y'); % Dodanie legendy w optymalnym miejscu

% 5. Dodanie siatki
grid on; % Włączenie siatki


%###########################################################################################
%Reakcja R_2
%###########################################################################################
% 2. Tworzenie wykresu
figure; % Otwiera nowe okno wykresu
plot(z_plot, R_2, '-g', 'LineWidth', 3); 
hold on; 
plot(z_plot, R_2x, '-r', 'LineWidth', 1); 
hold on; 
plot(z_plot, R_2y, '-b', 'LineWidth', 1);
hold off;

% 3. Dodanie opisu osi i tytułu
xlabel('szerokosc przedmiotu z [mm]'); % Opis osi X
ylabel('Siła reakcji R2 i jej składowe [N]'); % Opis osi Y
title('Siła reakcji R_2'); % Tytuł wykresu

% 4. Dodanie legendy
legend('R_2', 'R_2x', 'R_2y'); % Dodanie legendy w optymalnym miejscu

% 5. Dodanie siatki
grid on; % Włączenie siatki

%###########################################################################################
%Reakcja R_3
%###########################################################################################
% 2. Tworzenie wykresu
figure; % Otwiera nowe okno wykresu
plot(z_plot, R_3, '-g', 'LineWidth', 3); 
hold on; 
plot(z_plot, R_3x, '-r', 'LineWidth', 1); 
hold on; 
plot(z_plot, R_3y, '-b', 'LineWidth', 1);
hold off;

% 3. Dodanie opisu osi i tytułu
xlabel('szerokosc przedmiotu z [mm]'); % Opis osi X
ylabel('Siła reakcji R_3 i jej składowe [N]'); % Opis osi Y
title('Siła reakcji R_3'); % Tytuł wykresu

% 4. Dodanie legendy
legend('R_3', 'R_3x', 'R_3y'); % Dodanie legendy w optymalnym miejscu

% 5. Dodanie siatki
grid on; % Włączenie siatki

%###########################################################################################
%Siła do siłownika
%###########################################################################################
% 2. Tworzenie wykresu
figure; % Otwiera nowe okno wykresu
plot(z_plot, -S, '-g', 'LineWidth', 3);
hold off;

% 3. Dodanie opisu osi i tytułu
xlabel('szerokosc przedmiotu z [mm]'); % Opis osi X
ylabel('S [N]'); % Opis osi Y
title('Siła generowana na siłowniku'); % Tytuł wykresu

% 4. Dodanie legendy
%legend('R_3', 'R_3x', 'R_3y'); % Dodanie legendy w optymalnym miejscu

% 5. Dodanie siatki
grid on; % Włączenie siatki