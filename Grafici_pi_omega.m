clear 
close all
clc


r_mc=(1.02725 + 0.41090)/4;
r_hp=(0.60930 + 0.43919)/4;
omega_mc=5112.79;
ns_mc=8;
omega_hp=12000.08;
ns_hp=6;
cp=1005;
gamma=1.4;
alpha_1=.3; %scelta saggiamente 
cz=1.290596727248291e+02; %velocità assiale che si conserva in tutto il compressore.

P_fan=1.266484859530714; %ottenuta da rapporti_compressione di Lorenzo.
Tt_1=2.791438848106087e+02;  % ottenuta da rapporti_compressione di Lorenzo.




Pi=zeros(130,130);

%%
k = 1;
h = 1;
for i=1:100:13000
    omega_mc=(i-1)*(2*pi/60);
    k = 1;

    for j=1:100:13000
        omega_hp=(j-1)*(2*pi/60);       
        
        Pi_vect=[1]; %comincio sempre da 1, non dal rapporto di compress del fan
        Tt_vect=[Tt_1];
        n=1;      
        %compressore mc
        while n<=ns_mc
            
            Pi_mc_new=(1+.9*(omega_mc*r_mc*(omega_mc*r_mc-2*cz*tan(alpha_1)))/(cp*Tt_vect(end)))^(gamma/(gamma-1));
            Tt_new=Tt_vect(end)+omega_mc*r_mc*(omega_mc*r_mc-2*cz*tan(alpha_1))/cp; 
            n=n+1;
            Pi_vect=[Pi_vect Pi_mc_new*Pi_vect(end)];
            Tt_vect=[Tt_vect Tt_new];
        end
        %compressore hp
         n=1;
        while n<=ns_hp
            Pi_hp_new=(1+.9*(omega_hp*r_hp*(omega_hp*r_hp-2*cz*tan(alpha_1)))/(cp*Tt_vect(end)))^(gamma/(gamma-1));
             Tt_new=Tt_vect(end)+omega_hp*r_hp*(omega_hp*r_hp-2*cz*tan(alpha_1))/cp; 
            n=n+1;
            Pi_vect=[Pi_vect Pi_hp_new*Pi_vect(end)];
            Tt_vect=[Tt_vect Tt_new];
        end
        % lungo le colonne aumento la omeg_hp
        Pi(h,k)=Pi_vect(end);
        k = k+1;

    end
    %lungo le righe aumento la omega_mc
    h = h+1;
end
%%
omega_mc_vect=[1:100:13000];
omega_hp_vect=[1:100:13000];
wMC = 5112; %rpm
wHP = 12000.08; %rpm


figure(2)
surf(omega_mc_vect,omega_hp_vect,Pi)
title('rapporto di compressione a diverse velocità di rotazione');
xlabel('velocità di rotazioe MC [RPM]');
ylabel('velocità di rotazioe HP [RPM]');
zlabel('rapporto di compressione');
hold on
plot3(wMC,wHP,Pi(120,51),'w .','markersize',70)
hold off
grid on
        
        
    
    
    
    
    
    
    
    
    
    
    