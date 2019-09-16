
clear
clc
close all

%%

global low_index upp_index gamma_xy_fit w skip win_index L_fit Phi_fit

%% Da li se koriste rezultati dobijeni u FRESPID ili NAVFIT modulu?

uneto = 0; 

while(uneto==0)
    prethodni_modul = input('\n Prethodni modul je: 1 - FRESPID, 2 - COMPOSITE, 3 - MISOSA! \n');
    uneto = pravilan_unos(prethodni_modul,1,3);
end

%% 

% ako se koriste rezultati iz FRESPID modula...

if (prethodni_modul == 1)
    
    
    ime_fajla = input('\n Unesi ime .mat fajla (npr: "ime.mat"): \n');
    % ucitaj podatke
    load([ime_fajla '_frespid' '.mat'])
    
    k = 0;
    for i = 1:broj_razlichitih_prozora
 %% Prikaz procenjene uchestanosne karakteristike

    % ... vrshi se prikaz dobijenih rezultata za vishe prozora
    
        figure(12)
        subplot(3,1,1)
        semilogx(w,L_Amp(:,i),'Color',[0/255 k/255 255/255])
        hold on
        
            if (i==broj_razlichitih_prozora)

                i1 = find(w<=0.8*wmin,1,'last');
                i2 = find(w>=1.2*wmax,1,'first');
                Ylow = min(min( L_Amp(i1:i2,:) ));
                Yhigh = max(max( L_Amp(i1:i2,:) ));

                if (sign(Yhigh)==1)
                    Yhigh = 1.2*Yhigh;
                else
                    Yhigh = 0.8*Yhigh;
                end
                
                if (sign(Ylow)==1)
                    Ylow = 0.8*Ylow;
                else
                    Ylow = 1.2*Ylow;
                end
                
                axis([0.8*wmin 1.2*wmax Ylow Yhigh])
                grid on
                hold on

                xlabel('$\omega [\mathrm{rad/s}]$','Interpreter','Latex','FontSize',14);
                ylabel('$ L(\omega) [\mathrm{dB}]$','Interpreter','Latex','FontSize',14);
                legend({'$a = 0.25$', '$a = 0.3$', '$ a = 0.35$'}, 'Interpreter','Latex', 'Location', 'southwest')
            end

            subplot(3,1,2)
            semilogx(w,Phi(:,i),'Color',[0/255 k/255 255/255])
            hold on
    
            if (i==broj_razlichitih_prozora)

                i1 = find(w<=0.8*wmin,1,'last');
                i2 = find(w>=1.2*wmax,1,'first');
                Ylow = min(min( Phi(i1:i2,:) ));
                Yhigh = max(max( Phi(i1:i2,:) ));
                
                if (sign(Yhigh)==1)
                    Yhigh = 1.2*Yhigh;
                else
                    Yhigh = 0.8*Yhigh;
                end
                
                if (sign(Ylow)==1)
                    Ylow = 0.8*Ylow;
                else
                    Ylow = 1.2*Ylow;
                end
                
                axis([0.8*wmin 1.2*wmax Ylow Yhigh])
                grid on
                hold on

                xlabel('$\omega [\mathrm{rad/s}]$','Interpreter','Latex','FontSize',14);
                ylabel('$ \varphi(\omega) [^{\circ}] $','Interpreter','Latex','FontSize',14);
                legend({'$a = 0.25$', '$a = 0.3$', '$ a = 0.35$'}, 'Interpreter','Latex', 'Location', 'southwest')
            end

            subplot(3,1,3)
            semilogx(w,gamma_xy(:,i),'Color',[255/255 k/255 0])
            hold on
    
            if (i==broj_razlichitih_prozora)

                i1 = find(w<=0.8*wmin,1,'last');
                i2 = find(w>=1.2*wmax,1,'first');
                Ylow = min(min( gamma_xy(i1:i2,:) ));
                Yhigh = max(max( gamma_xy(i1:i2,:) ));
                
                if (sign(Yhigh)==1)
                    Yhigh = 1.2*Yhigh;
                else
                    Yhigh = 0.8*Yhigh;
                end
                
                if (sign(Ylow)==1)
                    Ylow = 0.8*Ylow;
                else
                    Ylow = 1.2*Ylow;
                end
                
                axis([0.8*wmin 1.2*wmax Ylow Yhigh])
                grid on

                xlabel('$\omega [\mathrm{rad/s}]$','Interpreter','Latex','FontSize',14);
                ylabel('$ \gamma_{xy}^2 $','Interpreter','Latex','FontSize',14);
                semilogx([0.8*wmin 1.2*wmin],[0.6 0.6],'k--')
                legend({'$a = 0.25$', '$a = 0.3$', '$ a = 0.35$'}, 'Interpreter','Latex', 'Location', 'southwest')
                set(gcf,'Position',[10 -200 1000 1200])
                
            end
            
        k = k + 255/broj_razlichitih_prozora;
        
    end
    
    % ...bira se zheljena procena od svih dostupnih koja se odredi kao
    % relevantna

    uneto = 0; 

    while(uneto==0)
        win_index = input('\n Unesi redni broj prozora koji uzimas kao relevantan! \n');
        uneto = pravilan_unos(win_index,1,broj_razlichitih_prozora);
    end

%% Za odabrani prozor (procenu), prikazhi sve!

    for i = win_index
 %% Prikaz procenjene uchestanosne karakteristike

        figure(13)
        subplot(3,1,1)
        semilogx(w,L_Amp(:,i),'b')
        hold on
    
        i1 = find(w<=0.8*wmin,1,'last');
        i2 = find(w>=1.2*wmax,1,'first');
        
        Ylow = min(min( L_Amp(i1:i2,win_index) ));
        Yhigh = max(max( L_Amp(i1:i2,win_index) ));

        if (sign(Yhigh)==1)
            Yhigh = 1.2*Yhigh;
        else
            Yhigh = 0.8*Yhigh;
        end
        
        if (sign(Ylow)==1)
            Ylow = 0.8*Ylow;
        else
            Ylow = 1.2*Ylow;
        end
        
        axis([0.8*wmin 1.2*wmax Ylow Yhigh])
        grid on
        hold on

        xlabel('$\omega [\mathrm{rad/s}]$','Interpreter','Latex','FontSize',14);
        ylabel('$ L(\omega) [\mathrm{dB}]$','Interpreter','Latex','FontSize',14);

        subplot(3,1,2)
        semilogx(w,Phi(:,i),'b')
        hold on
    
        i1 = find(w<=0.8*wmin,1,'last');
        i2 = find(w>=1.2*wmax,1,'first');

        Ylow = min(min( Phi(i1:i2,win_index) ));
        Yhigh = max(max( Phi(i1:i2,win_index) ));
    
        if (sign(Yhigh)==1)
           Yhigh = 1.2*Yhigh;
        else
           Yhigh = 0.8*Yhigh;
        end
        
        if (sign(Ylow)==1)
            Ylow = 0.8*Ylow;
        else
            Ylow = 1.2*Ylow;
        end
    
        axis([0.8*wmin 1.2*wmax Ylow Yhigh])
        grid on
        hold on

        xlabel('$\omega [\mathrm{rad/s}]$','Interpreter','Latex','FontSize',14);
        ylabel('$ \varphi(\omega) [^{\deg}] $','Interpreter','Latex','FontSize',14);

        subplot(3,1,3)
        semilogx(w,gamma_xy(:,i),'b')
        hold on

        i1 = find(w<=0.8*wmin,1,'last');
        i2 = find(w>=1.2*wmax,1,'first');

        Ylow = min(min( gamma_xy(i1:i2,win_index) ));
        Yhigh = max(max( gamma_xy(i1:i2,win_index) ));
    
        if (sign(Yhigh)==1)
            Yhigh = 1.2*Yhigh;
        else
            Yhigh = 0.8*Yhigh;
        end
        
        if (sign(Ylow)==1)
            Ylow = 0.8*Ylow;
        else
            Ylow = 1.2*Ylow;
        end
        axis([0.8*wmin 1.2*wmax Ylow Yhigh])
        grid on
    
        xlabel('$\omega [\mathrm{rad/s}]$','Interpreter','Latex','FontSize',14);
        ylabel('$ \gamma_xy^2 $','Interpreter','Latex','FontSize',14);
        semilogx([0.8*wmin 1.2*wmin],[0.6 0.6],'k--')
        set(gcf,'Position',[10 -200 1000 1200])

    end
 
    % ako je odabran izlaz iz modula COMPOSITE...
    
elseif (prethodni_modul==2)
    
    ime_fajla = input('Unesi ime .mat fajla (npr: "ime.mat"): \n');
    load([ime_fajla '_composite' '.mat'])
    
    win_index = 1;
    
    for i = win_index
 %% Prikaz procenjene uchestanosne karakteristike
        
 % ... prikazhi dobijene procene iz COMPOSITE modula
        figure(13)
        subplot(3,1,1)
        semilogx(w,L_Amp(:,i),'b')
        hold on
    
        i1 = 1;
        i2 = length(w);
        
        Ylow = min(min( L_Amp(i1:i2,win_index) ));
        Yhigh = max(max( L_Amp(i1:i2,win_index) ));

        if (sign(Yhigh)==1)
            Yhigh = 1.2*Yhigh;
        else
            Yhigh = 0.8*Yhigh;
        end
        
        if (sign(Ylow)==1)
            Ylow = 0.8*Ylow;
        else
            Ylow = 1.2*Ylow;
        end
        axis([0.8*wmin 1.2*wmax Ylow Yhigh])
        grid on
        hold on

        xlabel('$\omega [\mathrm{rad/s}]$','Interpreter','Latex','FontSize',14);
        ylabel('$ L(\omega) [\mathrm{dB}] $','Interpreter','Latex','FontSize',14);

        subplot(3,1,2)
        semilogx(w,Phi(:,i),'b')
        hold on
    
        i1 = 1;
        i2 = length(w);

        Ylow = min(min( Phi(i1:i2,win_index) ));
        Yhigh = max(max( Phi(i1:i2,win_index) ));
    
        if (sign(Yhigh)==1)
           Yhigh = 1.2*Yhigh;
        else
           Yhigh = 0.8*Yhigh;
        end
        
        if (sign(Ylow)==1)
            Ylow = 0.8*Ylow;
        else
            Ylow = 1.2*Ylow;
        end
    
        axis([0.8*wmin 1.2*wmax Ylow Yhigh])
        grid on
        hold on

        xlabel('$\omega [\mathrm{rad/s}]$','Interpreter','Latex','FontSize',14);
        ylabel('$ \varphi(\omega) [^{\deg}] $','Interpreter','Latex','FontSize',14);

        subplot(3,1,3)
        semilogx(w,gamma_xy(:,i),'b')
        hold on

        i1 = 1;
        i2 = length(w);

        Ylow = min(min( gamma_xy(i1:i2,win_index) ));
        Yhigh = max(max( gamma_xy(i1:i2,win_index) ));
    
        if (sign(Yhigh)==1)
            Yhigh = 1.2*Yhigh;
        else
            Yhigh = 0.8*Yhigh;
        end
        
        if (sign(Ylow)==1)
            Ylow = 0.8*Ylow;
        else
            Ylow = 1.2*Ylow;
        end
        
        axis([0.8*wmin 1.2*wmax Ylow Yhigh])
        grid on
    
        xlabel('$\omega [\mathrm{rad/s}]$','Interpreter','Latex','FontSize',14);
        ylabel('$ \gamma_xy^2 $','Interpreter','Latex','FontSize',14);
        semilogx([0.8*wmin 1.2*wmin],[0.6 0.6],'k--')
        set(gcf,'Position',[10 -200 1000 850])
        
    end
    
else
    
    ime_fajla = input('Unesi ime .mat fajla (npr: "ime.mat"): \n');
    load([ime_fajla '_misosa' '.mat'])
    
    
uneto = 0; 
[~,nu] = size(x_detrend_fil);

while(uneto==0)
    broj_ulazne = input('Unesi broj ulazne velichine po kojoj se trazhi PF: \n');
    uneto = pravilan_unos(broj_ulazne,1,nu);
end
    
   
    
    win_index = broj_ulazne;
 %% Prikaz procenjene uchestanosne karakteristike
        
 % ... prikazhi dobijene procene iz MISOSA modula
        figure(13)
        subplot(3,1,1)
        semilogx(w,L_Amp(:,win_index),'b')
        hold on
    
        i1 = 1;
        i2 = length(w);
        
        Ylow = min(min( L_Amp(i1:i2,win_index) ));
        Yhigh = max(max( L_Amp(i1:i2,win_index) ));

        if (sign(Yhigh)==1)
            Yhigh = 1.2*Yhigh;
        else
            Yhigh = 0.8*Yhigh;
        end
        
        if (sign(Ylow)==1)
            Ylow = 0.8*Ylow;
        else
            Ylow = 1.2*Ylow;
        end
        axis([0.8*wmin 1.2*wmax Ylow Yhigh])
        grid on
        hold on

        xlabel('$\omega [\mathrm{rad/s}]$','Interpreter','Latex','FontSize',14);
        ylabel('$ L(\omega) [\mathrm{dB}] $','Interpreter','Latex','FontSize',14);

        subplot(3,1,2)
        semilogx(w,Phi(:,win_index),'b')
        hold on
    
        i1 = 1;
        i2 = length(w);

        Ylow = min(min( Phi(i1:i2,win_index) ));
        Yhigh = max(max( Phi(i1:i2,win_index) ));
    
        if (sign(Yhigh)==1)
           Yhigh = 1.2*Yhigh;
        else
           Yhigh = 0.8*Yhigh;
        end
        
        if (sign(Ylow)==1)
            Ylow = 0.8*Ylow;
        else
            Ylow = 1.2*Ylow;
        end
    
        axis([0.8*wmin 1.2*wmax Ylow Yhigh])
        grid on
        hold on

        xlabel('$\omega [\mathrm{rad/s}]$','Interpreter','Latex','FontSize',14);
        ylabel('$ \varphi(\omega) [^{\deg}] $','Interpreter','Latex','FontSize',14);

        subplot(3,1,3)
        semilogx(w,gammay_ddot_q(:,1),'b')
        hold on

        i1 = 1;
        i2 = length(w);

        Ylow = min(min( gammay_ddot_q(i1:i2,1) ));
        Yhigh = max(max( gammay_ddot_q(i1:i2,1) ));
    
        if (sign(Yhigh)==1)
            Yhigh = 1.2*Yhigh;
        else
            Yhigh = 0.8*Yhigh;
        end
        
        if (sign(Ylow)==1)
            Ylow = 0.8*Ylow;
        else
            Ylow = 1.2*Ylow;
        end
        
        axis([0.8*wmin 1.2*wmax Ylow Yhigh])
        grid on
    
        xlabel('$\omega [\mathrm{rad/s}]$','Interpreter','Latex','FontSize',14);
        ylabel('$ \gamma_xy^2 $','Interpreter','Latex','FontSize',14);
        semilogx([0.8*wmin 1.2*wmax],[0.6 0.6],'k--')
        set(gcf,'Position',[10 -200 1000 850])
        
  
    
end

 %% Do koje uchestanosti se fituje i u koliko tachaka?

% Na osnovu poslednjeg prikaza uneti maksimalnu uchestanost do koje se
% fituje

uneto = 0; 

while(uneto==0)
    w_max_fit = input('\n Unesi maksimalnu uchestanost do koje se vrshi fitovanje! \n');
    uneto = pravilan_unos(w_max_fit,wmin,wmax);
end

% indeksi minimalne i maksimalne uchestanosti za fitovanje
low_index = i1;
upp_index = i2;

% broj tachaka u kriterijumu optimalnosti
uneto = 0; 
while(uneto==0)
broj_tacaka_za_fit = input('\n Unesi broj tachaka za fitovanje: \n');
uneto = pravilan_unos(broj_tacaka_za_fit,1,length(low_index:upp_index));
end

duz = length(low_index:1:upp_index);

skip = floor(duz/broj_tacaka_za_fit);

L_fit = L_Amp(low_index:skip:upp_index,win_index);
Phi_fit = Phi(low_index:skip:upp_index,win_index)*pi/180;

if (prethodni_modul == 1)
gamma_xy_fit = gamma_xy(low_index:skip:upp_index,1);
end

if (prethodni_modul == 3)
gamma_xy_fit = gammay_ddot_q(low_index:skip:upp_index,1);
end


%% Prikaz pochetnog reshenja

%x0 = [10^(-3.5/20);0.8];
x0 = [28;0.02;0.06;0.05]
%x0 = [10^(-4/20);1]
%x0 = [0.63;0.5]
%x0 = [25;3.5;1.05;0.1;1082];

figure(13)
W0 = pretpostavljena(x0)
[AMP,PHI] = bode(W0,w);
subplot(3,1,1)
semilogx(w,20*log10(AMP(1,:)),'k')
subplot(3,1,2)
semilogx(w,PHI(1,:),'k')

%%
xm = sign(x0).*(abs(x0)/20);
xM = sign(x0).*(20*abs(x0));

for i = 1:length(xm)
    if(xm(i) <0);
        a = xm(i);
        xm(i) = xM(i);
        xM(i) = a;
    end
end

%% Odredjivanje optimalnih parametara modela

opcije = optimset('Algorithm','interior-point','TolX',1e-6,'TolFun',1e-6,'MaxFunEvals',10000);

xopt = fmincon(@(x) opt_fun(x),x0,[],[],[],[],xm,xM,[[],[]],opcije);

J = opt_fun(xopt);


fprintf('Vrednost kriterijuma optimalnosti je: %f \n', J)

%% Prikaz optimalnog reshenja
% s = tf('s');
% figure(13)
% West1 = pretpostavljena(xopt)
% [AMP,PHI] = bode(West1,w);
% subplot(3,1,1)
% semilogx(w,20*log10(AMP(1,:)),'c')
% subplot(3,1,2)
% semilogx(w,PHI(1,:),'c')

%%
xm = sign(xopt).*(abs(xopt)/20);
xM = sign(xopt).*(20*abs(xopt));

for i = 1:length(xm)
    if(xm(i) <0);
        a = xm(i);
        xm(i) = xM(i);
        xM(i) = a;
    end
end

opcije = optimset('Algorithm','interior-point','TolX',1e-4,'TolFun',1e-4,'MaxFunEvals',10000);

xopt = fmincon(@(x) opt_fun(x),3*xopt,[],[],[],[],xm,xM,[[],[]],opcije);

J = opt_fun(xopt);

fprintf('Vrednost kriterijuma optimalnosti je: %f \n', J)

%%
%% Prikaz optimalnog reshenja

% figure(13)
% West2 = pretpostavljena(xopt)
% [AMP,PHI] = bode(West2,w);
% subplot(3,1,1)
% semilogx(w,20*log10(AMP(1,:)),'m')
% subplot(3,1,2)
% semilogx(w,PHI(1,:),'m')

%%
xm = sign(xopt).*(abs(xopt)/20);
xM = sign(xopt).*(20*abs(xopt));

for i = 1:length(xm)
    if(xm(i) <0);
        a = xm(i);
        xm(i) = xM(i);
        xM(i) = a;
    end
end

opcije = optimset('Algorithm','sqp','TolX',1e-6,'TolFun',1e-6,'MaxFunEvals',10000);

xopt = fmincon(@(x) opt_fun(x),0.5*xopt,[],[],[],[],xm,xM,[[],[]],opcije);

J = opt_fun(xopt);

fprintf('Vrednost kriterijuma optimalnosti je: %f \n', J)

%%
%% Prikaz optimalnog reshenja

figure(13)
West3 = pretpostavljena(xopt)
[AMP,PHI] = bode(West3,w);
subplot(3,1,1)
semilogx(w,20*log10(AMP(1,:)),'r')
legend({'$\hat{W}$', '$W_0$', '$W$'}, 'Interpreter','Latex')
subplot(3,1,2)
semilogx(w,PHI(1,:),'r')
legend({'$\hat{W}$', '$W_0$', '$W$'}, 'Interpreter','Latex', 'Location', 'southwest')
set(gcf,'Position',[10 -200 900 1200])

