clc; clear;
L_1=680*10^(-3);
L_2=0.8*L_1;
a_1=0.4*L_1;
a_2=0.3*L_2;
J_2=1.05;
J_1=J_2/0.8;
alfa=70*pi/180;
fi_bis=5.01;
Q_2=95;
Q_1=Q_2/0.9;
G=100;
a=270*10^(-3);
b=140*10^(-3);
h=200*10^(-3);
r_h=380*10^(-3);
g=9.81;
omega=[0 -1;
       1 0];
fi_stopnie=20:0.01:100;
n=length(fi_stopnie);
for i=1:1:n
    fi=fi_stopnie(i)*pi/180;
    r_c1=a_1*[cos(fi) sin(fi)]';
    r_s=Rot(fi)*[r_h h]';
    r_b=[-a b]';
    
    gamma=atan2((r_s(2)-r_b(2))/(norm(r_s-r_b)),(r_s(1)-r_b(1))/(norm(r_s-r_b)));
    gamma_stopnie=gamma*180/pi;
    r_c2=Rot(fi)*([L_1 0]'+Rot(-(pi/2+(pi/2-alfa)))*[a_2 0]');
    
    r_G=Rot(fi)*([L_1 0]'+Rot(-(pi/2+(pi/2-alfa)))*[L_2 0]');
    
    m_1=Q_1/g;
    m_2=Q_2/g;
    
    J_10=J_1+m_1*a_1*a_1;
    J_20=J_2+m_2*norm(r_c2)*norm(r_c2);
    J_G0=G*norm(r_G)*norm(r_G)/g;
    J_0=J_10+J_20+J_G0;
    Q_1v=[0 -Q_1]';
    Q_2v=[0 -Q_2]';
    Gv=[0 -G]';
    S(i)=(+J_0*fi_bis-[omega*r_c1]'*Q_1v-[omega*r_c2]'*Q_2v-[omega*r_G]'*Gv)/([omega*r_s]'*[cos(gamma) sin(gamma)]');
    S_x(i)=S(i)*cos(gamma);
    S_y(i)=S(i)*sin(gamma);

    eta=0.85;
    M=2.39;
    S_odM=-2000*pi*eta/10;

    R_x(i)=-S_odM*cos(gamma);
    R_y(i)=G+Q_1+Q_2-S_odM*sin(gamma);
    R(i)=(R_x(i)*R_x(i)+R_y(i)*R_y(i))^(0.5);


    %do śrub w łączniku prostopapadłymi i ramieniu
    r_Sm=30*sqrt(2)*10^(-3);
    beta(i)=90-fi;
    S_M=S_odM*sin(beta+gamma);
    S_1(i)=S_odM/4;
    S_2(i)=0.25*S_M(i)*h/r_Sm;

    S_w11(1:2,i)=S_1(i)*[cos(beta(i)+gamma) sin(beta(i)+gamma)]'+S_2(i)*[cos(5*pi/4+pi) sin(5*pi/4+pi)]';
    S_w21(1:2,i)=S_1(i)*[cos(beta(i)+gamma) sin(beta(i)+gamma)]'+S_2(i)*[cos(7*pi/4+pi) sin(7*pi/4+pi)]';
    S_w31(1:2,i)=S_1(i)*[cos(beta(i)+gamma) sin(beta(i)+gamma)]'+S_2(i)*[cos(1*pi/4+pi) sin(1*pi/4+pi)]';
    S_w41(1:2,i)=S_1(i)*[cos(beta(i)+gamma) sin(beta(i)+gamma)]'+S_2(i)*[cos(3*pi/4+pi) sin(3*pi/4+pi)]';

    S_w11_norm(i)=norm(S_w11(1:2,i));
    S_w21_norm(i)=norm(S_w21(1:2,i));
    S_w31_norm(i)=norm(S_w31(1:2,i));
    S_w41_norm(i)=norm(S_w41(1:2,i));
end

%###########################################################################################
%Siła S
%###########################################################################################
% 2. Tworzenie wykresu
figure; % Otwiera nowe okno wykresu
plot(fi_stopnie, S, '-g', 'LineWidth', 3);
hold on; 
plot(fi_stopnie, S_x, '-r', 'LineWidth', 1);
hold on;
plot(fi_stopnie, S_y, '-b', 'LineWidth', 1);
hold off;

% 3. Dodanie opisu osi i tytułu
xlabel('fi [°]'); % Opis osi X
ylabel('S [N]'); % Opis osi Y
title('Siła S - ruch przeciwnie do ruchu wskazówek zegara'); % Tytuł wykresu
legend('S', 'S_x', 'S_y')
% 4. Dodanie legendy
%legend('R_3', 'R_3x', 'R_3y'); % Dodanie legendy w optymalnym miejscu

% 5. Dodanie siatki
grid on; % Włączenie siatki

%###########################################################################################
%Siła R
%###########################################################################################
% 2. Tworzenie wykresu
figure; % Otwiera nowe okno wykresu
plot(fi_stopnie, R, '-g', 'LineWidth', 3);

hold off;

% 3. Dodanie opisu osi i tytułu
xlabel('fi [°]'); % Opis osi X
ylabel('R [N]'); % Opis osi Y
title('Siła reakcji w podstawie - siła S działa w stronę do silnika'); % Tytuł wykresu

% 4. Dodanie legendy
%legend('R_3', 'R_3x', 'R_3y'); % Dodanie legendy w optymalnym miejscu

% 5. Dodanie siatki
grid on; % Włączenie siatki


%#####################################################################################
%Siły w śrubach w łączniku prostopadłym
%#####################################################################################
% 2. Tworzenie wykresu
figure; % Otwiera nowe okno wykresu
plot(fi_stopnie, S_w11_norm, '-g', 'LineWidth', 1);
hold on; 
plot(fi_stopnie, S_w21_norm, '-r', 'LineWidth', 1);
hold on;
plot(fi_stopnie, S_w31_norm, '-y', 'LineWidth', 1);
hold on;
plot(fi_stopnie, S_w41_norm, '-b', 'LineWidth', 1);
hold off;

% 3. Dodanie opisu osi i tytułu
xlabel('fi [°]'); % Opis osi X
ylabel('S_w [N]'); % Opis osi Y
title('Siła wypadkowe w śrubach w łączniku przy ramieniu'); % Tytuł wykresu
legend('S_{w1}', 'S_{w2}', 'S_{w3}', 'S_{w4}');
% 4. Dodanie legendy
%legend('R_3', 'R_3x', 'R_3y'); % Dodanie legendy w optymalnym miejscu

% 5. Dodanie siatki
grid on; % Włączenie siatki



for i=1:1:n
    fi=fi_stopnie(i)*pi/180;
    r_c1=a_1*[cos(fi) sin(fi)]';
    r_s=Rot(fi)*[r_h h]';
    r_b=[-a b]';
    
    gamma=atan2((r_s(2)-r_b(2))/(norm(r_s-r_b)),(r_s(1)-r_b(1))/(norm(r_s-r_b)));
    gamma_stopnie=gamma*180/pi;
    r_c2=Rot(fi)*([L_1 0]'+Rot(-(pi/2+(pi/2-alfa)))*[a_2 0]');
    
    r_G=Rot(fi)*([L_1 0]'+Rot(-(pi/2+(pi/2-alfa)))*[L_2 0]');
    
    m_1=Q_1/g;
    m_2=Q_2/g;
    
    J_10=J_1+m_1*a_1*a_1;
    J_20=J_2+m_2*norm(r_c2)*norm(r_c2);
    J_G0=G*norm(r_G)*norm(r_G)/g;
    J_0=J_10+J_20+J_G0;
    Q_1v=[0 -Q_1]';
    Q_2v=[0 -Q_2]';
    Gv=[0 -G]';
    S(i)=(-J_0*fi_bis-[omega*r_c1]'*Q_1v-[omega*r_c2]'*Q_2v-[omega*r_G]'*Gv)/([omega*r_s]'*[cos(gamma) sin(gamma)]');
    S_x(i)=S(i)*cos(gamma);
    S_y(i)=S(i)*sin(gamma);
    
    eta=0.85;
    M=2.39;
    S_odM=2000*pi*eta/10;

    R_x(i)=-S_odM*cos(gamma);
    R_y(i)=G+Q_1+Q_2-S_odM*sin(gamma);
    R(i)=(R_x(i)*R_x(i)+R_y(i)*R_y(i))^(0.5);

    %do śrub w łączniku prostopapadłymi i ramieniu
    r_Sm=30*sqrt(2)*10^(-3);
    beta(i)=90-fi;
    S_M=S_odM*sin(beta+gamma);
    S_1(i)=S_odM/4;
    S_2(i)=0.25*S_M(i)*h/r_Sm;

    S_w11(1:2,i)=S_1(i)*[cos(beta(i)+gamma) sin(beta(i)+gamma)]'+S_2(i)*[cos(5*pi/4+pi) sin(5*pi/4+pi)]';
    S_w21(1:2,i)=S_1(i)*[cos(beta(i)+gamma) sin(beta(i)+gamma)]'+S_2(i)*[cos(7*pi/4+pi) sin(7*pi/4+pi)]';
    S_w31(1:2,i)=S_1(i)*[cos(beta(i)+gamma) sin(beta(i)+gamma)]'+S_2(i)*[cos(1*pi/4+pi) sin(1*pi/4+pi)]';
    S_w41(1:2,i)=S_1(i)*[cos(beta(i)+gamma) sin(beta(i)+gamma)]'+S_2(i)*[cos(3*pi/4+pi) sin(3*pi/4+pi)]';

    S_w11_norm(i)=norm(S_w11(1:2,i));
    S_w21_norm(i)=norm(S_w21(1:2,i));
    S_w31_norm(i)=norm(S_w31(1:2,i));
    S_w41_norm(i)=norm(S_w41(1:2,i));
end

%###########################################################################################
%Siła S
%###########################################################################################
% 2. Tworzenie wykresu
figure; % Otwiera nowe okno wykresu
plot(fi_stopnie, S, '-g', 'LineWidth', 3);
hold on; 
plot(fi_stopnie, S_x, '-r', 'LineWidth', 1);
hold on;
plot(fi_stopnie, S_y, '-b', 'LineWidth', 1);
hold off;

% 3. Dodanie opisu osi i tytułu
xlabel('fi [°]'); % Opis osi X
ylabel('S [N]'); % Opis osi Y
title('Siła S - ruch zgodnie ze wskazówkami zegara'); % Tytuł wykresu
legend('S', 'S_x', 'S_y')
% 4. Dodanie legendy
%legend('R_3', 'R_3x', 'R_3y'); % Dodanie legendy w optymalnym miejscu

% 5. Dodanie siatki
grid on; % Włączenie siatki

%###########################################################################################
%Siła R
%###########################################################################################
% 2. Tworzenie wykresu
figure; % Otwiera nowe okno wykresu
plot(fi_stopnie, R, '-g', 'LineWidth', 3);

hold off;

% 3. Dodanie opisu osi i tytułu
xlabel('fi [°]'); % Opis osi X
ylabel('R [N]'); % Opis osi Y
title('Siła reakcji w podstawie - siła od śruby działa w stronę przeciwną do silnika'); % Tytuł wykresu

% 4. Dodanie legendy
%legend('R_3', 'R_3x', 'R_3y'); % Dodanie legendy w optymalnym miejscu

% 5. Dodanie siatki
grid on; % Włączenie siatki


%#####################################################################################
%Siły w śrubach w łączniku prostopadłym
%#####################################################################################
% 2. Tworzenie wykresu
figure; % Otwiera nowe okno wykresu
plot(fi_stopnie, S_w11_norm, '-g', 'LineWidth', 1);
hold on; 
plot(fi_stopnie, S_w21_norm, '-r', 'LineWidth', 1);
hold on;
plot(fi_stopnie, S_w31_norm, '-y', 'LineWidth', 1);
hold on;
plot(fi_stopnie, S_w41_norm, '-b', 'LineWidth', 1);
hold off;

% 3. Dodanie opisu osi i tytułu
xlabel('fi [°]'); % Opis osi X
ylabel('S_w [N]'); % Opis osi Y
title('Siła wypadkowe w śrubach w łączniku przy ramieniu'); % Tytuł wykresu
legend('S_{w1}', 'S_{w2}', 'S_{w3}', 'S_{w4}');
% 4. Dodanie legendy
%legend('R_3', 'R_3x', 'R_3y'); % Dodanie legendy w optymalnym miejscu

% 5. Dodanie siatki
grid on; % Włączenie siatki


%##########################################################
%SPRZĘGŁO
%##########################################################

R_e=235*10^6;%Stal S235JR
v=0.3;
E=2.1*10^11;
d=5*10^(-3);
r=23.75*10^(-3);
M=8.36;
k=1.5;


F=M/(4*r);
n=400;
x = linspace(0, 18, n);
for i=1:n
    if x(i)<3
        Mg_x(i)=0;
        T_x(i)=0;
    end
    if x(i)>=3 && x(i)<9
        Mg_x(i)=-0.5*F*(x(i)-3)*10^(-3);
        T_x(i)=-0.5*F;
    end
    if x(i)>=9 && x(i)<15
        Mg_x(i)=-0.5*F*(x(i)-3)*10^(-3)+F*(x(i)-9)*10^(-3);
        T_x(i)=-0.5*F+F;   
    end
    if x(i)>=15 && x(i)<=18
        Mg_x(i)=-0.5*F*(x(i)-3)*10^(-3)+F*(x(i)-9)*10^(-3)-0.5*F*(x(i)-15)*10^(-3);
        T_x(i)=-0.5*F+F-0.5*F;   
    end
end

    % 2. Tworzenie wykresu
figure; % Otwiera nowe okno wykresu
plot(x, Mg_x, '-g', 'LineWidth', 3);
hold off;

% 3. Dodanie opisu osi i tytułu
xlabel('x [mm]'); % Opis osi X
ylabel('Mg [Nm]'); % Opis osi Y
title('Moment gnący w kołku sprzęgła'); % Tytuł wykresu

% 5. Dodanie siatki
grid on; % Włączenie siatki

    % 2. Tworzenie wykresu
figure; % Otwiera nowe okno wykresu
plot(x, T_x, '-g', 'LineWidth', 3);
hold off;

% 3. Dodanie opisu osi i tytułu
xlabel('x [mm]'); % Opis osi X
ylabel('T [N]'); % Opis osi Y
title('Siła tnąca w kołku sprzęgła'); % Tytuł wykresu

% 5. Dodanie siatki
grid on; % Włączenie siatki
M=-0.263008;
T=-44.0;

sigma=(32*M/(pi*(d^3)));
tau=16*T/(3*pi*d*d);
sigma_red=sqrt(sigma^2 + 3* tau^2);
sigma_max=k*sigma_red;

disp('#############################################');
disp('#############################################');
disp('Sprzęgło - kołki - naprężenia');
disp('#############################################');
disp('F=')
disp(F)
disp('sigma=')
disp(sigma)
disp('tau=')
disp(tau)
disp('sigma_red=')
disp(sigma_red)
disp('sigma_max=')
disp(sigma_max)
if((sigma_max)<R_e)
    disp('Kołki wytrzymają');
else 
    disp('Kołki nie wytrzymają');
end
disp('#############################################');
disp('#############################################');
disp(' ');

%#######################################################
%SPRZĘGLO Kołek - ścinanie
%#######################################################
m=1;
tau=F/(0.25*pi*d*d*m);
disp('#############################################');
disp('#############################################');
disp('Sprzęgło - kołki - ścinanie');
disp('#############################################');
disp('sigma=')
disp('tau')
disp('k*tau=')
disp(k*tau)
if((k*tau)<R_e)
    disp('Kołki wytrzymają ścinanie');
else 
    disp('Kołki nie wytrzymają ścinanie');
end
disp('#############################################');
disp('#############################################');
disp(' ');



%#######################################################
%SPRZĘGLO ŚRUBY
%#######################################################
S=133.8; 
Q_x=S/4;
E_k=2.1*10^11; %Eelementy sprzęgła wykonane ze stali 42CrMo4
E_kp=2.1*10^11; %Podkładka: Stal C45
E_s=2.1*10^11;
l_s=18*10^(-3);
l_k=6*10^(-3);
d=5*10^(-3);
s=8*10^(-3);
g=6*10^(-3);
d_3=4.019*10^(-3);
k_r=210*10^6; %stal A2 AISI 304 i jej granica sprężystości


A_s=(0.8*d)*(0.8*d)*pi/4;
A_k=(s+g)*(s+g)*pi/4;
tan_alfa=A_s*E_s/l_s
tan_beta=1/((l_k/A_k)*((2/E_k)+(1/E_kp)))
Q_w=Q_x*(1-(tan_alfa/tan_beta));
F_kr=pi*d_3*d_3*k_r/4;

disp('#############################################');
disp('Sprzeglo - sruby');
disp('#############################################');
disp('naciag wstepny:');
disp(Q_w);
disp('Q_x:');
disp(Q_x);
disp('k*Q_x:');
disp(k*Q_x);
disp('F_kr:');
disp(F_kr);
if(k*Q_x<F_kr)
    disp('Śruby M5 wytrzymają');
else 
    disp('Śruby M5 nie wytrzymają');
end

%#######################################################
%Łącznik prostopadły - śruby
%#######################################################


R_e=500*10^6; %Stal 25CrMo4 i jej granica sprężystości
Fx=1298.46; %odczytane z wykresu dla śruby 2
F=Fx;
n=400;
x = linspace(0, 104.5, n);
for i=1:n
    if x(i)<8
        Mg_x(i)=0;
        T_x(i)=0;
    end
    if x(i)>=8 && x(i)<51
        Mg_x(i)=-0.5*F*(x(i)-8)*10^(-3);
        T_x(i)=-0.5*F;
    end
    if x(i)>=51 && x(i)<94
        Mg_x(i)=-0.5*F*(x(i)-8)*10^(-3)+F*(x(i)-51)*10^(-3);
        T_x(i)=-0.5*F+F;   
    end
    if x(i)>=94 && x(i)<=104.5
        Mg_x(i)=-0.5*F*(x(i)-8)*10^(-3)+F*(x(i)-51)*10^(-3)-0.5*F*(x(i)-94)*10^(-3);
        T_x(i)=-0.5*F+F-0.5*F;   
    end
end


    % 2. Tworzenie wykresu
figure; % Otwiera nowe okno wykresu
plot(x, Mg_x, '-g', 'LineWidth', 3);
hold off;

% 3. Dodanie opisu osi i tytułu
xlabel('x [mm]'); % Opis osi X
ylabel('Mg [Nm]'); % Opis osi Y
title('Moment gnący w śrubie nr. 2'); % Tytuł wykresu

% 5. Dodanie siatki
grid on; % Włączenie siatki

    % 2. Tworzenie wykresu
figure; % Otwiera nowe okno wykresu
plot(x, T_x, '-g', 'LineWidth', 3);
hold off;

% 3. Dodanie opisu osi i tytułu
xlabel('x [mm]'); % Opis osi X
ylabel('T [N]'); % Opis osi Y
title('Siła tnąca w w śrubie nr. 2'); % Tytuł wykresu
grid on;
M=-27.8705; %spisane z wykresu
T=-649.23; %spisane z wykresu
d_3=11.546*10^(-3);
k=1.5;

sigma=(32*M/(pi*(d_3^3)));
tau=16*T/(3*pi*d_3*d_3);
sigma_red=sqrt(sigma^2 + 3* tau^2);
sigma_max=k*sigma_red;

disp('#############################################');
disp('#############################################');
disp('Łącznik prostopadły - śruby');
disp('#############################################');
disp('ZGINANIE')
disp('sigma=')
disp(sigma)
disp('tau=')
disp(tau)
disp('sigma_red=')
disp(sigma_red)
disp('sigma_max=')
disp(sigma_max)
if((sigma_max)<R_e)
    disp('śruby wytrzymają Zginanie');
else 
    disp('śruby nie wytrzymają Zginania');
end

disp('#############################################');
disp('#############################################');
disp(' ');



%#######################################################
%#######################################################


%PODSTAWA

%#######################################################
%#######################################################

%Momenty gnące w wale
n=1000;
R_e=700*10^6; %Stal 42CrMo4
R_max=834.5;
d=17*10^(-3);
x = linspace(0, 141.5, n);


for i=1:n
    if x(i)<16.1
        Mg_x(i)=0;
        T_x(i)=0;
    end
    if x(i)>=16.1 && x(i)<70.75
        Mg_x(i)=-0.5*R_max*(x(i)-16.1)*10^(-3);
        T_x(i)=-0.5*R_max;
    end
    if x(i)>=70.75 && x(i)<125.4
        Mg_x(i)=-0.5*R_max*(x(i)-16.1)*10^(-3)+R_max*(x(i)-70.75)*10^(-3);
        T_x(i)=0.5*R_max;  
    end
    if x(i)>=125.4 && x(i)<=141.5
        Mg_x(i)=0;
        T_x(i)=0;   
    end
end


    % 2. Tworzenie wykresu
figure; % Otwiera nowe okno wykresu
plot(x, Mg_x, '-g', 'LineWidth', 3);
hold off;

% 3. Dodanie opisu osi i tytułu
xlabel('x [mm]'); % Opis osi X
ylabel('Mg [Nm]'); % Opis osi Y
title('Moment gnący w wale'); % Tytuł wykresu

% 5. Dodanie siatki
grid on; % Włączenie siatki

    % 2. Tworzenie wykresu
figure; % Otwiera nowe okno wykresu
plot(x, T_x, '-g', 'LineWidth', 3);
hold off;

% 3. Dodanie opisu osi i tytułu
xlabel('x [mm]'); % Opis osi X
ylabel('T [N]'); % Opis osi Y
title('Siła tnąca w wale'); % Tytuł wykresu
grid on;

M=-22.7732;
T=-417.25;
sigma=(32*M/(pi*(d^3)));
tau=16*T/(3*pi*d*d);
sigma_red=sqrt(sigma^2 + 3* tau^2);
sigma_max=k*sigma_red;

disp('#############################################');
disp('#############################################');
disp('PODSTAWA - WAŁ - naprężenia');
disp('#############################################');
disp('sigma=')
disp(sigma)
disp('tau=')
disp(tau)
disp('sigma_red=')
disp(sigma_red)
disp('sigma_max=')
disp(sigma_max)
if((sigma_max)<R_e)
    disp('Wał wytrzyma');
else 
    disp('Wał wytrzyma');
end
disp('#############################################');
disp('#############################################');
disp(' ');




%######################################################################
%Łożysko SKF32303
%######################################################################
disp('#############################################');
disp('#############################################');
disp('Łożysko SKF32303');
disp('#############################################');
C_0=33.5*10^3;
Y_0=1.1;
P_p=0.5*R_max
%tutaj kąt obliczono uprzednio na podstawie modelu CAD, ale wynik pokrył
%się z tym, co umieścił producent na stronie łożyska takż nie sprawdzano
%już tego nie zmieniano
gamma2=9.6822*pi/180;
delta2=1.0758*pi/180;
9.6822+1.0758
P_w=P_p*tan(gamma2+delta2/2)

stosunekPw_C0=P_w/C_0;

stosunekPw_Pp=P_w/P_p

P_0=0.5*P_p+Y_0*P_w

if P_0<P_p
        disp('P_0=P_p');
        P_0=P_p
else
    disp('P_0=0.5*P_p+Y_0*P_w');
    P_0=0.5*P_p+Y_0*P_w
end

if P_0<C_0
        disp('Łożysko SKF32303 wytrzyma');

else
    disp('Łożysko SKF32303 nie wytrzyma');
end

disp(' ');
disp('Łożysko SKF32303 ale siła P_p jest przemnożona przez k');
disp('#############################################');
k=1.5;

C_0=33.5*10^3;
Y_0=1.1;
P_p=k*0.5*R_max
gamma2=9.6822*pi/180;
delta2=1.0758*pi/180;

P_w=P_p*tan(gamma2+delta2/2)

stosunekPw_C0=P_w/C_0;

stosunekPw_Pp=P_w/P_p

P_0=0.5*P_p+Y_0*P_w

if P_0<P_p
        disp('P_0=P_p');
        P_0=P_p
else
    disp('P_0=0.5*P_p+Y_0*P_w');
    P_0=0.5*P_p+Y_0*P_w
end

if P_0<C_0
        disp('Łożysko SKF32303 wytrzyma');

else
    disp('Łożysko SKF32303 nie wytrzyma');
end

disp('#############################################');
disp('#############################################');
disp('WYBOCZENIE');
disp('#############################################');
disp('#############################################');
f_k=0.0625
d_r=16.6
L_t=586.3
F_k=4.072*10^5*((f_k*d_r^4)/(L_t^2))
F_p=0.5*F_k
k=1.5
S_Tp=1276.4 
S_max=k*S_Tp

