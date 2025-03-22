clc; clear;
L=40;
z_plot=30:0.01:52.5; %%%%Zakres wymiarów chwytanego przedmiotu dzielony na 2
P_max=344;
a=15;
b=22;
d=79;
c=10;
m=8.0;

omega =[0 -1; 
        1 0];

vector_Pmax=[-P_max/2 0]';
vector_c=[0 10]';
vector_r3=[0 -(b-c)]';

z_plot=2*z_plot;
i=1;
for z=30.0:0.01:52.5
    
    rn=52.5+m;

    alfa(i) = asin((d-(z+c))/L);
    beta(i) = -asin(((z+c)-(a+b))/L);
    
    vector_rn=[rn 0]';
    
     b_vector=[0; 0; 0; -P_max/(2*cos(alfa(i)))];
     
     A=[cos(alfa(i)) cos(beta(i)) cos(beta(i)) 0;
        sin(alfa(i)) sin(beta(i)) sin(beta(i)) 1;
        (omega*vector_c)'*[cos(alfa(i)) sin(alfa(i))]' (omega*vector_r3)'*[cos(beta(i)) sin(beta(i))]' (omega*vector_rn)'*[cos(beta(i)) sin(beta(i))]' [omega*vector_rn]' *[0 1]';
        1 0 0 0];
     
     
     q=A\b_vector;
%###########################################################################################
%DANE DO WYKRESÓW
%###########################################################################################

%siły reakcji

     R_1(i)=q(1);
     R_2(i)=q(2);
     R_3(i)=q(3);
     N(i)=q(4);
     
     R_1x(i)=R_1(i)*cos(alfa(i));
     R_1y(i)=R_1(i)*sin(alfa(i));

     R_2x(i)=R_2(i)*cos(beta(i));
     R_2y(i)=R_2(i)*sin(beta(i));
     
     R_3x(i)=R_3(i)*cos(beta(i));
     R_3y(i)=R_3(i)*sin(beta(i));
    
    %siła na siłownik
    S(i)=2*R_1(i)*cos(alfa(i));


    i=i+1;
end




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
xlabel('szerokość przedmiotu z [mm]'); % Opis osi X
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
xlabel('szerokość przedmiotu z [mm]'); % Opis osi X
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
xlabel('szerokość przedmiotu z [mm]'); % Opis osi X
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
plot(z_plot, 2*S, '-g', 'LineWidth', 3);
hold off;

% 3. Dodanie opisu osi i tytułu
xlabel('szerokość przedmiotu z [mm]'); % Opis osi X
ylabel('S [N]'); % Opis osi Y
title('Siła generowana na siłowniku'); % Tytuł wykresu

% 4. Dodanie legendy
%legend('R_3', 'R_3x', 'R_3y'); % Dodanie legendy w optymalnym miejscu

% 5. Dodanie siatki
grid on; % Włączenie siatki


%###########################################################################################
%Siła nacisku
%###########################################################################################
% 2. Tworzenie wykresu
figure; % Otwiera nowe okno wykresu
plot(z_plot, N, '-g', 'LineWidth', 3);
hold off;

% 3. Dodanie opisu osi i tytułu
xlabel('szerokość przedmiotu z [mm]'); % Opis osi X
ylabel('N [N]'); % Opis osi Y
title('Siła nacisku na przedmiot'); % Tytuł wykresu

% 4. Dodanie legendy
%legend('R_3', 'R_3x', 'R_3y'); % Dodanie legendy w optymalnym miejscu

% 5. Dodanie siatki
grid on; % Włączenie siatki

%###########################################################################################
%MOMENT GNĄCY I SIŁA TNĄCA
%###########################################################################################
%materiał:  Stal S1100QL
display('Stal S1100QL')
E=2.1*10^(11)
v=0.3
d=6*10^(-3);
Z_g=1100*10^(6) %Wytrzymałość na zginanie w zakresie 1000-1200MPa ale przyjmuje najmniej korzystny scenariusz
k=1.5;
z_max=3*10^(-3);
I_y=(pi*d^(4))/64;


%Sworzeń 3
Rx_3=R_3(1)/2;
n=400;
x_3 = linspace(0, 58, n);

for i=1:n
    if x_3(i)<6
        Mg_x3(i)=0;
        T_x3(i)=0;
    end
    if x_3(i)>=6 && x_3(i)<17.5
        Mg_x3(i)=-Rx_3*(x_3(i)-6)*10^(-3);
        T_x3(i)=-Rx_3;
    end
    if x_3(i)>=17.5 && x_3(i)<40.5
        Mg_x3(i)=-Rx_3*(x_3(i)-6)*10^(-3)+Rx_3*(x_3(i)-17.5)*10^(-3);
        T_x3(i)=0;
    end
    if x_3(i)>=40.5 && x_3(i)<52
        Mg_x3(i)=-Rx_3*(x_3(i)-6)*10^(-3)+Rx_3*(x_3(i)-17.5)*10^(-3)+Rx_3*(x_3(i)-40.5)*10^(-3);
        T_x3(i)= Rx_3;
    end
    if x_3(i)>=52 && x_3(i)<=58
        Mg_x3(i)=-Rx_3*(x_3(i)-6)*10^(-3)+Rx_3*(x_3(i)-17.5)*10^(-3)+Rx_3*(x_3(i)-40.5)*10^(-3)-Rx_3*(x_3(i)-52)*10^(-3);
        T_x3(i)=0;
    end  
end

%n
%length(Mg_x3)
%length(T_x3)
    % 2. Tworzenie wykresu
figure; % Otwiera nowe okno wykresu
plot(x_3, Mg_x3, '-g', 'LineWidth', 3);
hold off;

% 3. Dodanie opisu osi i tytułu
xlabel('x [mm]'); % Opis osi X
ylabel('Mg [Nm]'); % Opis osi Y
title('Moment gnący w sworzniu nr. 3'); % Tytuł wykresu

% 5. Dodanie siatki
grid on; % Włączenie siatki

    % 2. Tworzenie wykresu
figure; % Otwiera nowe okno wykresu
plot(x_3, T_x3, '-g', 'LineWidth', 3);
hold off;

% 3. Dodanie opisu osi i tytułu
xlabel('x [mm]'); % Opis osi X
ylabel('T [N]'); % Opis osi Y
title('Siła tnąca w sworzniu nr. 3'); % Tytuł wykresu

% 5. Dodanie siatki
grid on; % Włączenie siatki

M=-17.4993;
T=-1521.68;

sigma=(32*M/(pi*(d^3)));
tau=16*T/(3*pi*d*d);

sigma_red=sqrt(0.5*sigma^2 + 3* tau^2);

sigma_max=k*sigma_red;

disp('#############################################');
disp('#############################################');
disp('Sworzeń 3');
disp('#############################################');
disp('sigma=')
disp(sigma)
disp('tau=')
disp(tau)
disp('sigma_red=')
disp(sigma_red)
disp('sigma_max=')
disp(sigma_max)
if((sigma_max)<Z_g)
    disp('Sworzeń 3 wytrzyma');
else 
    disp('Sworzeń 3 nie wytrzyma');
end
disp('#############################################');
disp('#############################################');
disp(' ');

%Sworzeń 2
R_x2=R_1(1)/2;
x_2 = linspace(0, 30, n);

for i=1:n
    if x_2(i)<5
        Mg_x2(i)=0;
        T_x2(i)=0;
    end
    if x_2(i)>=5 && x_2(i)<15
        Mg_x2(i)=R_x2*(x_2(i)-5)*10^(-3);
        T_x2(i)=R_x2;
    end
    if x_2(i)>=15 && x_2(i)<25
        Mg_x2(i)=R_x2*(x_2(i)-5)*10^(-3)-2*R_x2*(x_2(i)-15)*10^(-3);
        T_x2(i)=-R_x2;
    end
    if x_2(i)>=25 && x_2(i)<=30
        Mg_x2(i)=R_x2*(x_2(i)-5)*10^(-3)-2*R_x2*(x_2(i)-15)*10^(-3)+R_x2*(x_2(i)-25)*10^(-3);
        T_x2(i)=0;
    end

end

    % 2. Tworzenie wykresu
figure; % Otwiera nowe okno wykresu
plot(x_2, Mg_x2, '-g', 'LineWidth', 3);
hold off;

% 3. Dodanie opisu osi i tytułu
xlabel('x [mm]'); % Opis osi X
ylabel('Mg [Nm]'); % Opis osi Y
title('Moment gnący w sworzniu nr. 2'); % Tytuł wykresu

% 5. Dodanie siatki
grid on; % Włączenie siatki

    % 2. Tworzenie wykresu
figure; % Otwiera nowe okno wykresu
plot(x_2, T_x2, '-g', 'LineWidth', 3);
hold off;

% 3. Dodanie opisu osi i tytułu
xlabel('x [mm]'); % Opis osi X
ylabel('T [N]'); % Opis osi Y
title('Siła tnąca w sworzniu nr. 2'); % Tytuł wykresu

% 5. Dodanie siatki
grid on; % Włączenie siatki

M=-3.85575;
T=-387.03;

sigma=(32*M/(pi*(d^3)));
tau=16*T/(3*pi*d*d);

sigma_red=sqrt(0.5*sigma^2 + 3* tau^2);

sigma_max=k*sigma_red;

disp('#############################################');
disp('#############################################');
disp('Sworzeń 2');
disp('#############################################');
disp('sigma=')
disp(sigma)
disp('tau=')
disp(tau)
disp('sigma_red=')
disp(sigma_red)
disp('sigma_max=')
disp(sigma_max)
if((sigma_max)<Z_g)
    disp('Sworzeń 2 wytrzyma');
else 
    disp('Sworzeń 2 nie wytrzyma');
end
disp('#############################################');
disp('#############################################');
disp(' ');

%Sworzeń 4
R_x4=R_3(1)/2;
x_4 = linspace(0, 24, n);
R_R2=-R_x4*(12-6)/(23-12);
R_R1=R_x4-R_R2;
for i=1:n
    if x_4(i)<6
        Mg_x4(i)=0;
        T_x4(i)=0;
    end
    if x_4(i)>=6 && x_4(i)<12
        Mg_x4(i)=-R_x4*(x_4(i)-6)*10^(-3);
        T_x4(i)=-R_x4;
    end
    if x_4(i)>=12 && x_4(i)<=23
        Mg_x4(i)=-R_x4*(x_4(i)-6)*10^(-3)+R_R1*(x_4(i)-12)*10^(-3);
        T_x4(i)=-R_x4+R_R1;
    end
    if x_4(i)>=23 && x_4(i)<=24
        Mg_x4(i)=-R_x4*(x_4(i)-6)*10^(-3)+R_R1*(x_4(i)-12)*10^(-3)+R_R2*(x_4(i)-23)*10^(-3);
        T_x4(i)=0;
    end

end

    % 2. Tworzenie wykresu
figure; % Otwiera nowe okno wykresu
plot(x_4, Mg_x4, '-g', 'LineWidth', 3);
hold off;

% 3. Dodanie opisu osi i tytułu
xlabel('x [mm]'); % Opis osi X
ylabel('Mg [Nm]'); % Opis osi Y
title('Moment gnący w sworzniu nr. 4'); % Tytuł wykresu

% 5. Dodanie siatki
grid on; % Włączenie siatki

    % 2. Tworzenie wykresu
figure; % Otwiera nowe okno wykresu
plot(x_4, T_x4, '-g', 'LineWidth', 3);
hold off;

% 3. Dodanie opisu osi i tytułu
xlabel('x [mm]'); % Opis osi X
ylabel('T [N]'); % Opis osi Y
title('Siła tnąca w sworzniu nr. 4'); % Tytuł wykresu

% 5. Dodanie siatki
grid on; % Włączenie siatki

M=-9.10509;
T=-1521.68;

sigma=(32*M/(pi*(d^3)));
tau=16*T/(3*pi*d*d);

sigma_red=sqrt(0.5*sigma^2 + 3* tau^2);

sigma_max=k*sigma_red;

disp('#############################################');
disp('#############################################');
disp('Sworzeń 4');
disp('#############################################');
disp('sigma=')
disp(sigma)
disp('tau=')
disp(tau)
disp('sigma_red=')
disp(sigma_red)
disp('sigma_max=')
disp(sigma_max)
if((sigma_max)<Z_g)
    disp('Sworzeń 4 wytrzyma');
else 
    disp('Sworzeń 4 nie wytrzyma');
end
disp('#############################################');
disp('#############################################');
disp(' ');


%Sworzeń 5
R_x5=R_2(1)/2;
x_5 = linspace(0, 24, n);
R_R2=-R_x5*(12-6)/(23-12);
R_R1=R_x5-R_R2;
for i=1:n
    if x_5(i)<6
        Mg_x5(i)=0;
        T_x5(i)=0;
    end
    if x_5(i)>=6 && x_5(i)<12
        Mg_x5(i)=-R_x5*(x_5(i)-6)*10^(-3);
        T_x5(i)=-R_x5;
    end
    if x_5(i)>=12 && x_5(i)<=23
        Mg_x5(i)=-R_x5*(x_5(i)-6)*10^(-3)+R_R1*(x_5(i)-12)*10^(-3);
        T_x5(i)=-R_x5+R_R1;
    end
    if x_5(i)>=23 && x_5(i)<=24
        Mg_x5(i)=-R_x5*(x_5(i)-6)*10^(-3)+R_R1*(x_5(i)-12)*10^(-3)+R_R2*(x_5(i)-23)*10^(-3);
        T_x5(i)=0;
    end

end

    % 2. Tworzenie wykresu
figure; % Otwiera nowe okno wykresu
plot(x_5, Mg_x5, '-g', 'LineWidth', 3);
hold off;

% 3. Dodanie opisu osi i tytułu
xlabel('x [mm]'); % Opis osi X
ylabel('Mg [Nm]'); % Opis osi Y
title('Moment gnący w sworzniu nr. 5'); % Tytuł wykresu

% 5. Dodanie siatki
grid on; % Włączenie siatki

    % 2. Tworzenie wykresu
figure; % Otwiera nowe okno wykresu
plot(x_5, T_x5, '-g', 'LineWidth', 3);
hold off;

% 3. Dodanie opisu osi i tytułu
xlabel('x [mm]'); % Opis osi X
ylabel('T [N]'); % Opis osi Y
title('Siła tnąca w sworzniu nr. 5'); % Tytuł wykresu

% 5. Dodanie siatki
grid on; % Włączenie siatki

M=8.58905;
T=1435.43;

sigma=(32*M/(pi*(d^3)));
tau=16*T/(3*pi*d*d);

sigma_red=sqrt(0.5*sigma^2 + 3* tau^2);

sigma_max=k*sigma_red;
disp('#############################################');
disp('#############################################');
disp('Sworzeń 5');
disp('#############################################');
disp('sigma=')
disp(sigma)
disp('tau=')
disp(tau)
disp('sigma_red=')
disp(sigma_red)
disp('sigma_max=')
disp(sigma_max)
if((sigma_max)<Z_g)
    disp('Sworzeń 5 wytrzyma');
else 
    disp('Sworzeń 5 nie wytrzyma');
end
disp('#############################################');
disp('#############################################');
disp(' ');



%Sworzeń 1

%oś x

Rx_1=R_1x(1);
Rx_2=0.5*R_2x(1);
R_R2=(0.5*Rx_1*(29-12)-Rx_2*(12-6))/(23-12);
R_R1=Rx_2+0.5*Rx_1-R_R2;
n=1000;
z_1 = linspace(0, 58, n);

for i=1:n
    if z_1(i)<6
        Mgx_1(i)=0;
        Tx_1(i)=0;
    end
    if z_1(i)>=6 && z_1(i)<12;
        Mgx_1(i)=-Rx_2*(z_1(i)-6)*10^(-3);
        Tx_1(i)=-Rx_2;
    end
    if z_1(i)>=12 && z_1(i)<23;
        Mgx_1(i)=-Rx_2*(z_1(i)-6)*10^(-3)+R_R1*(z_1(i)-12)*10^(-3);
        Tx_1(i)=-Rx_2+R_R1;
    end
    if z_1(i)>=23 && z_1(i)<29
        Mgx_1(i)=-Rx_2*(z_1(i)-6)*10^(-3)+R_R1*(z_1(i)-12)*10^(-3)+R_R2*(z_1(i)-23)*10^(-3);
        Tx_1(i)=-Rx_2+R_R1+R_R2;
    end
    if z_1(i)>=29 && z_1(i)<35
        Mgx_1(i)=-Rx_2*(z_1(i)-6)*10^(-3)+R_R1*(z_1(i)-12)*10^(-3)+R_R2*(z_1(i)-23)*10^(-3)-Rx_1*(z_1(i)-29)*10^(-3);
        Tx_1(i)=-Rx_2+R_R1+R_R2-Rx_1;
    end
    if z_1(i)>=35 && z_1(i)<46
        Mgx_1(i)=-Rx_2*(z_1(i)-6)*10^(-3)+R_R1*(z_1(i)-12)*10^(-3)+R_R2*(z_1(i)-23)*10^(-3)-Rx_1*(z_1(i)-29)*10^(-3)+R_R2*(z_1(i)-35)*10^(-3);
        Tx_1(i)=-Rx_2+R_R1+R_R2-Rx_1+R_R2;
    end
    if z_1(i)>=46 && z_1(i)<=52
        Mgx_1(i)=-Rx_2*(z_1(i)-6)*10^(-3)+R_R1*(z_1(i)-12)*10^(-3)+R_R2*(z_1(i)-23)*10^(-3)-Rx_1*(z_1(i)-29)*10^(-3)+R_R2*(z_1(i)-35)*10^(-3)+R_R1*(z_1(i)-46)*10^(-3);
        Tx_1(i)=-Rx_2+R_R1+R_R2-Rx_1+R_R2+R_R1;
    end  
    if z_1(i)>=52 && z_1(i)<=58
        Mgx_1(i)=-Rx_2*(z_1(i)-6)*10^(-3)+R_R1*(z_1(i)-12)*10^(-3)+R_R2*(z_1(i)-23)*10^(-3)-Rx_1*(z_1(i)-29)*10^(-3)+R_R2*(z_1(i)-35)*10^(-3)+R_R1*(z_1(i)-46)*10^(-3)-Rx_2*(z_1(i)-52)*10^(-3);
        Tx_1(i)=-Rx_2+R_R1+R_R2-Rx_1+R_R2+R_R1-Rx_2;
    end  
end


%oś y
Ry_1=R_1y(1);
Ry_2=0.5*R_2y(1);
R_R2=(0.5*Ry_1*(29-12)-Ry_2*(12-6))/(23-12);
R_R1=Ry_2+0.5*Ry_1-R_R2;

for i=1:n
    if z_1(i)<6
        Mgy_1(i)=0;
        Ty_1(i)=0;
    end
    if z_1(i)>=6 && z_1(i)<12;
        Mgy_1(i)=-Ry_2*(z_1(i)-6)*10^(-3);
        Ty_1(i)=-Ry_2;
    end
    if z_1(i)>=12 && z_1(i)<23;
        Mgy_1(i)=-Ry_2*(z_1(i)-6)*10^(-3)+R_R1*(z_1(i)-12)*10^(-3);
        Ty_1(i)=-Ry_2+R_R1;
    end
    if z_1(i)>=23 && z_1(i)<29
        Mgy_1(i)=-Ry_2*(z_1(i)-6)*10^(-3)+R_R1*(z_1(i)-12)*10^(-3)+R_R2*(z_1(i)-23)*10^(-3);
        Ty_1(i)=-Ry_2+R_R1+R_R2;
    end
    if z_1(i)>=29 && z_1(i)<35
        Mgy_1(i)=-Ry_2*(z_1(i)-6)*10^(-3)+R_R1*(z_1(i)-12)*10^(-3)+R_R2*(z_1(i)-23)*10^(-3)-Ry_1*(z_1(i)-29)*10^(-3);
        Ty_1(i)=-Ry_2+R_R1+R_R2-Ry_1;
    end
    if z_1(i)>=35 && z_1(i)<46
        Mgy_1(i)=-Ry_2*(z_1(i)-6)*10^(-3)+R_R1*(z_1(i)-12)*10^(-3)+R_R2*(z_1(i)-23)*10^(-3)-Ry_1*(z_1(i)-29)*10^(-3)+R_R2*(z_1(i)-35)*10^(-3);
        Ty_1(i)=-Ry_2+R_R1+R_R2-Ry_1+R_R2;
    end
    if z_1(i)>=46 && z_1(i)<=52
        Mgy_1(i)=-Ry_2*(z_1(i)-6)*10^(-3)+R_R1*(z_1(i)-12)*10^(-3)+R_R2*(z_1(i)-23)*10^(-3)-Ry_1*(z_1(i)-29)*10^(-3)+R_R2*(z_1(i)-35)*10^(-3)+R_R1*(z_1(i)-46)*10^(-3);
        Ty_1(i)=-Ry_2+R_R1+R_R2-Ry_1+R_R2+R_R1;
    end  
    if z_1(i)>=52 && z_1(i)<=58
        Mgy_1(i)=-Ry_2*(z_1(i)-6)*10^(-3)+R_R1*(z_1(i)-12)*10^(-3)+R_R2*(z_1(i)-23)*10^(-3)-Ry_1*(z_1(i)-29)*10^(-3)+R_R2*(z_1(i)-35)*10^(-3)+R_R1*(z_1(i)-46)*10^(-3)-Ry_2*(z_1(i)-52)*10^(-3);
        Ty_1(i)=-Ry_2+R_R1+R_R2-Ry_1+R_R2+R_R1-Ry_2;
    end  
end


% 2. Tworzenie wykresu
figure; % Otwiera nowe okno wykresu
plot(z_1, Mgx_1, '-r', 'LineWidth', 1); 
hold off;

% 3. Dodanie opisu osi i tytułu
xlabel('z [mm]'); % Opis osi X
ylabel('Momenty gnące na poszczególnych osiach [Nm]'); % Opis osi Y
title('Wykres momentu gnącego na osi OX w sworzni nr. 1'); % Tytuł wykresu

% 4. Dodanie legendy
%legend('M_gx1', 'M_gy1'); % Dodanie legendy w optymalnym miejscu

% 5. Dodanie siatki
grid on; % Włączenie siatki

% 2. Tworzenie wykresu
figure; % Otwiera nowe okno wykresu
plot(z_1, Mgy_1, '-b', 'LineWidth', 1);
hold off;

% 3. Dodanie opisu osi i tytułu
xlabel('z [mm]'); % Opis osi X
ylabel('Mg  [Nm]'); % Opis osi Y
title('Wykres momentu gnącego na osi OY w sworzni nr. 1'); % Tytuł wykresu

% 4. Dodanie legendy
%legend('M_gx1', 'M_gy1'); % Dodanie legendy w optymalnym miejscu

% 5. Dodanie siatki
grid on; % Włączenie siatki

% 2. Tworzenie wykresu
figure; % Otwiera nowe okno wykresu
plot(z_1, Tx_1, '-r', 'LineWidth', 1); 
hold off;

% 3. Dodanie opisu osi i tytułu
xlabel('z [mm]'); % Opis osi X
ylabel('T [N]'); % Opis osi Y
title('Wykres siły tnącej na osi OX w sworzni nr. 1'); % Tytuł wykresu

% 4. Dodanie legendy
%legend('T_x1', 'T_y1'); % Dodanie legendy w optymalnym miejscu

% 5. Dodanie siatki
grid on; % Włączenie siatki

% 2. Tworzenie wykresu
figure; % Otwiera nowe okno wykresu
plot(z_1, Ty_1, '-b', 'LineWidth', 1);
hold off;

% 3. Dodanie opisu osi i tytułu
xlabel('z [mm]'); % Opis osi X
ylabel('T [N]'); % Opis osi Y
title('Wykres siły tnącej na osi OY w sworzni nr. 1'); % Tytuł wykresu

% 4. Dodanie legendy
%legend('T_x1', 'T_y1'); % Dodanie legendy w optymalnym miejscu

% 5. Dodanie siatki
grid on; % Włączenie siatki


%Wartości spisane ręcznie z wykresu w najbardziej wytężonym punkcie (17.5mm lub 40.5mm)
%Zaokrąglono wyniki do góry, aby naprężenia wyszły bardziej niekorzystne
Mx=8.57511;
My=-0.641178;
Tx=1431.39;
Ty=-107.657;

M=sqrt(Mx^2 + My^2);
T=sqrt(Tx^2 + Ty^2);

sigma=(32*M/(pi*(d^3)));
tau=16*T/(3*pi*d*d);

sigma_red=sqrt(0.5*sigma^2 + 3* tau^2);

sigma_max=k*sigma_red;
disp('#############################################');
disp('#############################################');
disp('Sworzeń 1');
disp('#############################################');
disp('Mg_1max=')
disp(M)
disp('T_1max')
disp(T)
disp('sigma=')
disp(sigma)
disp('tau=')
disp(tau)
disp('sigma_red=')
disp(sigma_red)
disp('sigma_max=')
disp(sigma_max)
if((sigma_max)<Z_g)
    disp('Sworzeń 1 wytrzyma');
else 
    disp('Sworzeń 1 nie wytrzyma');
end
disp('#############################################');
disp('#############################################');
disp(' ');

%############################################################################################
%Wytrzymałość Łożysk
%############################################################################################
disp(' ');
disp(' ');
disp(' ');
disp(' ');
disp(' ');
disp('#############################################');
disp('#############################################');
disp('Wytrzymałość łożysk');
disp('#############################################');
disp('#############################################');
disp(' ');
disp('Podane przez producenta maksymalne naprężenia powierzchniowe');
loz_sigma_max=80*10^(6)
disp('#############################################');
%łożysko GSM-0608-03
b1_3=3*10^(-3);
d1_3=6*10^(-3);
d2_3=8*10^(-3);
A_3=d2_3*b1_3

%łożysko GFM-0608-05
d1_5=6*10^(-3);
d2_5=8*10^(-3);
d3_5=12*10^(-3);
b1_5=5*10^(-3);
b2_5=1*10^(-3);
A_5=(b1_5-b2_5)*d2_5


%łącznik 1-2
disp('łącznik 1-2');
F=0.5*abs(R_1(1))
sigma=k*F/A_5
if sigma<loz_sigma_max
    disp("łożyska GFM-0608-05 wytrzymają");
end

disp(" ");
%łącznik 3-4
disp("####################################");
disp('łącznik 3-4');
F=0.25*abs(R_3(1))
sigma=k*F/A_5
if sigma<loz_sigma_max
    disp("łożyska GFM-0608-05 wytrzymają");
end

disp(" ");
%łącznik 1-5
disp("####################################");
disp('łącznik 1-5');
F=0.25*abs(R_2(1))
sigma=k*F/A_5
if sigma<loz_sigma_max
    disp("łożyska GFM-0608-05 wytrzymają");
end

disp(" ");
%rama
disp("####################################");
disp('rama');
F=0.25*abs(R_1(1))
sigma_5=k*F/A_5
if sigma_5<loz_sigma_max
    disp("łożyska GFM-0608-05 wytrzymają");
end
sigma_3=k*F/A_3
if sigma_3<loz_sigma_max
    disp("łożyska GSM-0608-03 wytrzymają");
end


disp(" ");
%palec sworzeń 1
disp("####################################");
disp('palec sworzeń 1');
F=0.25*norm(R_1(1)*[cos(alfa(1)) sin(alfa(i))]' + R_2(1)*[cos(beta(1)) sin(beta(i))]')

sigma_3=k*F/A_3
if sigma_3<loz_sigma_max
    disp("łożyska GSM-0608-03 wytrzymają");
end

disp(" ");
%palec sworzeń 3
disp("####################################");
disp('palec sworzeń 3');
F=0.25*abs(R_3(1))

sigma_5=k*F/A_5
if sigma_5<loz_sigma_max
    disp("łożyska GFM-0608-05 wytrzymają");
end
sigma_3=k*F/A_3
if sigma_3<loz_sigma_max
    disp("łożyska GSM-0608-03 wytrzymają");
end

disp(" ");
%ceownik sworzeń 4
disp("####################################");
disp('ceownik sworzeń 4');
F=0.25*abs(R_3(1))

sigma_5=k*F/A_5
if sigma_5<loz_sigma_max
    disp("łożyska GFM-0608-05 wytrzymają");
end
sigma_3=k*F/A_3
if sigma_3<loz_sigma_max
    disp("łożyska GSM-0608-03 wytrzymają");
end

disp(" ");
%ceownik sworzeń 5
disp("####################################");
disp('ceownik sworzeń 5');
F=0.25*abs(R_2(1))

sigma_5=k*F/A_5
if sigma_5<loz_sigma_max
    disp("łożyska GFM-0608-05 wytrzymają");
end
sigma_3=k*F/A_3
if sigma_3<loz_sigma_max
    disp("łożyska GSM-0608-03 wytrzymają");
end
