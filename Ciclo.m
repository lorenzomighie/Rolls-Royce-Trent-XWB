%Progetto di Propulsione, turbofan a flussi separati
clear
clc
%% Aria Standard

h=9000; %m

theta_0 = 288.15; %K
p_0 = 101325; % Pa
rho_0 = 1.225; %Kg/m^3
R = 287.05;
gamma=1.4;
a = sqrt(R*theta_0*gamma);

theta = theta_0-h*.0065; %K
p = p_0*(theta/theta_0)^(9.81/(.0065*R));
rho = rho_0*(theta/theta_0)^(-(1+9.81/(-.0065*R)));


%% Condizioni a monte   0

M_0=0.7; % scelto da noi
T_0=theta;
P_0=p;
rho_0=rho;

U=M_0*sqrt(gamma*R*T_0);
Pt_0=P_0*(1+((gamma-1)/2)*M_0^2)^(gamma/(gamma-1));
Tt_0=T_0*(1+((gamma-1)/2)*M_0^2);

%% Dati da algoritmo sugli stadi
Beta_fan = 1.266484859530714;
Beta_c = 47.115750605092920;
eta_fan = 0.896591465796934;
eta_c = 0.836861796623104;

%% Dati Motore Rolls Royce Trent 
pr_f=1.2;
pr_i=1.25; % 8 stadi
pr_h=1.4; % 6 stadi
D_fan=3.17944; % m
D_shaft=0.89978;
A_fan= pi*(D_fan/2)^2-pi*(D_shaft/2)^2;
BPR = 9.3; %di design, da file
f = 0.03;       
mm = 1442.923;

% Imponiamo Mach all'ingresso pari a .65 e all'ingresso del fan .4
M_fan=.54;
M_inlet=.65;
% la formula sotto vale per flusso isoentropico, poi però calcolo la caduta
% di pressione totale, con eps_presa, quindi dopo la considero non iso, ma
% solo adiabatica.

A_inlet=(M_fan/M_inlet)*A_fan*((1+((gamma-1)/2)*M_inlet^2)/(1+((gamma-1)/2)*M_fan^2))^((gamma+1)/(2*(gamma-1)));
D_inlet=sqrt(4*A_inlet/pi);

% Analogamente trovo la sezione del tubo di cattura, a mach .8, cioè di
% volo
% con _inlet indico sulla sezione di ingresso della presa

A_0=(M_inlet/M_0)*A_inlet*((1+((gamma-1)/2)*M_0^2)/(1+((gamma-1)/2)*M_inlet^2))^((gamma+1)/(2*(gamma-1)));
% con Bernulli ricavo la pressione statica all' inlet, compressione
% esterna, la totale rimane costante.

Tt_inlet=Tt_0;
T_inlet=Tt_inlet/(1+((gamma-1)/2)*M_inlet^2);
v_inlet=M_inlet*sqrt(gamma*R*T_inlet);
P_inlet=P_0+.5*rho_0*(U^2-v_inlet^2); % graficoooooooooooooooooooooo



%% Presa d'aria    1
% Supponiamo che la presa sia adiabatica, quindi la temperatura totale
% si conserva
% il rapporto di compressione si ottiene dal grafico della Sapienza, per
% mach di volo 0.8. Il rendimento adiabatico lo aggiungo ma in realtà non
% serve

Tt_1=Tt_0;
eps_presa=.96;
Pt_1=Pt_0* eps_presa;
eta_presa_adiabatico= (eps_presa^((gamma-1)/gamma)*(1+.5*(gamma-1)*M_0^2)-1)/(.5*(gamma-1)*M_0^2);

%la geometria della presa va da A_inlet a (A_fan+A_shaft)

%% presa d'aria 2

%continuità --> rapporto tra le aree
V0_inlet = a * M_inlet;
V1_inlet = a * M_fan;
A_inlet_inizio_r_Ainlet_fine = V1_inlet/V0_inlet;
A_inlet_exit = A_inlet/A_inlet_inizio_r_Ainlet_fine;
%teorema di Bernoulli (prima forma) HP di non viscosità --> ipotesi di
%Ptot = costante (serve solo per il CP)
P1_statica = (0.5*rho_0*V0_inlet^2 + p_0) - (rho_0*0.5*V1_inlet^2);

C_ppress = (P1_statica - p_0)/(0.5*rho_0*V0_inlet^2);
% evitare separazione 
%criterio Cp < 0.6

%angoli tra 5 e 7 gradi
% criterio grafico --> dimensiono 
theta_deg = 15;      %valore a STIMATO PER AVERE L GIUSTO O QUASI

theta_rad = theta_deg * pi / 180;
H = sqrt(A_inlet*4/pi); % altezza presa d'aria in igresso  --> fa venire
%un H troppo piccolo
B = sqrt(A_inlet_exit*4/pi); % altezza presa d'aria in uscita --> fa
%venire un B troppo piccolo

L = ((B-H)/2)/(tan(theta_rad/2)); % lunghezza del diffusore

% interpolazione grafico --> solo del transitorio --> evito anche la
% transizione, non solo la separazione così

Aspect_ratio =    [ 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 8, 9, 10]; 
%in gradi
Angolo_apertura = [ 15.5, 14, 13.2, 12.5, 11.9, 11, 10.4, 10, 9.7, 9, 8.2, 8];
    
xx = [3:0.1:10];
%figure
%plot (xx,spline(Aspect_ratio,Angolo_apertura, xx),'b');
%hold on
p = polyfit(Aspect_ratio,Angolo_apertura,5);
%plot (xx,polyval(p,xx),'r');
%grid on

%controllo sulla transizione
if (theta_deg < polyval(p,L/H))
    ok = 1;
else 
    ok = 0; 
end

%% Fan 2
% Beta_fan, il rapporto di compressione, è da ricavare o da cercare un valore tipico
% Qui i flussi non sono ancora separati
% Il rendimento adiabatico del fan è anche quello o da ricavare o da
% cercare

Pt_2=Pt_1*Beta_fan;
Tt_2_iso=Tt_1*Beta_fan^((gamma-1)/gamma);
Tt_2=Tt_1+(1/(eta_fan))*(Tt_2_iso-Tt_1);



%% Flusso aria fredda 5_cold
% Il BPR è da decidere o cercare
% Pi_fan indica le perdite lungo il condotto del fluido freddo.
Perdite_condotto=.975; % dal documento di Lore
Pt_5_cold= Perdite_condotto* Pt_2;
Tt_5_cold= Tt_2; % queste sono le condizioni all'uscita

%% Compressore 3
% Beta_c va calcolato dallo studio del comrpessore

Pt_3= Beta_c * Pt_2
Tt_3_iso = Tt_2 * Beta_c ^((gamma-1)/gamma);
Tt_3=Tt_2+(1/eta_c)*(Tt_3_iso-Tt_2);

Tt_compressore_medio=451.8243; %compressore MC


%% Camera di combustione  4
% f è da scegliere, così come Pi_b;
% gamma_gc e cp_gc sono anch'essi da determinare: CEA
% il rendimento delle camere di combustoine è di solito .99

%aggiungere diffusore per rallentare --> anche no
eta_b=0.9995; % valore preso da tabella

cp=R*gamma/(gamma-1);
Q_r=42e6 % J/kg
Pi_b=0.960;
cp_gc=1150;  %valori tipici dda 1145 a 1155

Tt_4=(1/((1+f)*cp_gc))*(cp*Tt_3+f*eta_b*Q_r);
Pt_4=Pt_3*Pi_b;

%% Turbina alta 5_high
% fa lavorare il compressore hp . Servono i rendimenti meccanici di
% turbina alta e compressore, eta_m_c e eta_m_th, e il rendimento
% adiabatico della turbina. eta_th
eta_m_c=.99;
% SERVE DIVIDERE MC DA HP
eta_m_th=.9;
gamma_gc=1.33;

eta_th = 0.902;
Tt_5_high= Tt_4-(cp/((1+f)*eta_m_c*eta_m_th*cp_gc))*( Tt_3 - Tt_compressore_medio);

Tt_5_high_iso = Tt_4-(1/eta_th)*(Tt_4-Tt_5_high);
Pt_5_high=Pt_4*(Tt_5_high_iso/Tt_4)^(gamma_gc/(gamma_gc-1));
%% Turbina meida 5_mid
% considero gli stessi rendimenti meccanici, nel file sono uguali
eta_tm = 0.9058;
Tt_5_mid=Tt_5_high-(cp/((1+f)*eta_m_c*eta_m_th*cp_gc))*( Tt_compressore_medio - Tt_2);

Tt_5_mid_iso = Tt_5_high-(1/eta_tm)*(Tt_5_high-Tt_5_mid);
Pt_5_mid=Pt_5_high*(Tt_5_mid_iso/Tt_5_high)^(gamma_gc/(gamma_gc-1));
%% Turbina bassa 5_low
% fa lavorare il fan: servono i rendimenti meccanici del fan e della
% turbina bassa, eta_m_f e eta_m_tl, e il rendimento adiabatico della
% turbina bassa eta_tl
eta_tl = 0.9189;
eta_m_f=.99;
eta_m_tl=.9;

Tt_5_low=Tt_5_mid-((1+BPR)*cp*(Tt_2-Tt_1))/((1+f)*eta_m_f*eta_m_tl*cp_gc); %sul fan ho sia aria calda che fredda

Tt_5_low_iso = Tt_5_mid-(1/eta_tl)*(Tt_5_mid-Tt_5_low);
Pt_5_low=Pt_5_mid*(Tt_5_low_iso/Tt_5_mid)^(gamma_gc/(gamma_gc-1));

%% Ugello freddo 6_cold
% Ipotesi di ugello ideale: si conservano sia la temperatura totale che la
% pressione totale. DA CONFERMARE.
% Ugello adattato, dunque la pressione statica è uguale a quella esterna

Pt_6_cold=Pt_5_cold;
Tt_6_cold=Tt_5_cold;
P_6_cold=P_0;

M_6_cold=sqrt( (2/(gamma-1))*((Pt_6_cold/P_6_cold)^((gamma-1)/gamma)-1));
T_6_cold=Tt_6_cold/(1+((gamma-1)/2)*M_6_cold^2);
a_6_cold=sqrt(gamma*R*T_6_cold);
v_6_cold=M_6_cold*a_6_cold;

%% Ugello caldo 6_hot

Pt_6_hot=Pt_5_low;
Tt_6_hot=Tt_5_low;
P_6_hot=P_0;

R_gc = cp_gc *(gamma_gc-1)/gamma_gc;
v_6_hot=sqrt( (2*gamma_gc/(gamma_gc-1))*R_gc*Tt_6_hot*(1-(P_6_hot/Pt_6_hot)^((gamma_gc-1)/gamma_gc)));
M_6_hot=sqrt( (2/(gamma_gc-1))*((Pt_6_hot/P_6_hot)^((gamma_gc-1)/gamma_gc)-1));

T_6_hot=Tt_6_hot/(1+((gamma_gc-1)/2)*M_6_hot^2);
a_6_hot=sqrt(gamma_gc*R*T_6_hot);


%% Spinta ecc
% la portata in ingresso deriva dai paramtri della presa d'aria e dalle
% condizioni di volo
% Le portate massiche sono mm_hot, mm_cold e mm_f, fuel, sono da definire
%definiti f e BPR --> basta dare una massa d'aria e si trova tutto!!

% in funzione di mm_hot
% BPR = mm_cold/mm_hot    e     f=mm_f/mm_hot
% T è la spinta
% mm è la portata aria in ingresso
% nota mm totale


mm_hot = mm/(BPR + 1);
mm_cold = mm_hot*BPR;
mm_f = f*mm_hot;

T= mm_hot*(1+f)*v_6_hot - mm_hot*U   +   mm_cold*(v_6_cold-U);

T = mm_hot*((1+f)*v_6_hot - (1+BPR)*U + BPR*v_6_cold)



TSFC=mm_f/T;

%% rendimento propulsivo
% effetto utile diviso effetto utile più l'energia cinetica persa
eta_prop = T*U/(T*U + 0.5*mm_hot*(v_6_hot-U)^2 + 0.5*mm_cold*(v_6_cold-U)^2);


%% rendimento termico

%potenza utile più energia cinetica del getto fratto potere calorifer
eta_termico = (T*U + 0.5*mm_hot*(v_6_hot-U)^2 + 0.5*mm_cold*(v_6_cold-U)^2)/(mm_f*Q_r);

%% rendimento globale

eta_globale = eta_termico*eta_prop;


%% formula dossi

eta_tdossih = ((1 + f)*v_6_hot^2 + BPR*v_6_cold^2 - (1+BPR)*U^2)/(2*f*Q_r);









