
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                        %
% ESTIMACIJA UCHESTANOSNE KARAKTERISTIKE JEDNOSTRUKO-PRENOSNOG SISTEMA   %
%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

clear
clc
close all
 
%% Broj nezavisnih snimaka

uneto = 0; 

while(uneto==0)
    % broj nezavisnih snimaka ponashanja sistema u okolini nominalne tachke
    N = input('\n Unesi broj nezavisnih snimaka: \n');
    uneto = pravilan_unos(N,1);
end

%% Spajanje ulaza i izlaza u vektore kolone

% vektor sa duzhinama trajanja pojedinih snimaka 
Trec = zeros(1,N);
% vektor sa periodama odabiranja u pojedinachnim snimcima
dT = zeros(1,N);
% vektor spojenih promena upravljanja
u_con = [];
% vektor spojenih promena nivoa 2 (izlazne velichine)
y_con = [];

%%

% za svaku od snimljenih promena...
for i = 1:N
    
    % reci kako se zove fajl sa podacima
    % (npr. 'snimak1.mat')
    ime_fajla = input('\n Unesi ime .mat fajla (npr: "ime.mat"): \n');
    
    % ucitaj podatke
    load(ime_fajla)
   
    % ovde voditi racuna o imenima izlaza i ulaza
    % (u originalnom snimku ulaz je u, a izlaz je y)
   
    % ako nije kolona, pretvori u kolonu
    [~,n] = size(u);
    
    if (n~=1)
        u = u';
    end
    
    % kako izlaz, tako i ulaz
    [~,n] = size(y);
    
    if (n~=1)
        y = y';
    end
   
    % pravi spojeni zapis spajanjem pojedinacnih 
    u_con = [u_con; u];
    y_con = [y_con; y];
   
    % odredi periodu odabiranja u pojedinacnim zapisima
    dT(1,i) = min(diff(time));
    
    % odredi duzhinu trajanja pojedinacnih zapisa
    Trec(1,i) = time(end);

end

% proveri da li su periode odabiranja iste za sve zapise
for i = 1:N
    
    if ( dT(i) ~= mean(dT) )
        error('\n Nisu iste ucestanosti odabiranja u snimanjima! \n')
    end
    
end

% odredi ucestanost odabiranja [Hz]
f_sampling  = 1/dT(1);

% i istu uchestanost u [rad/s]
w_sampling = 2*pi*f_sampling;

% ukupna duzhina trajana spojenih zapisa (TFinal)
TF = sum(Trec)+(N-1)*dT(1);

% vektor vremena koji sluzhi za prikaz spojenih zapisa
time_con = 0:dT(1):TF;

% obrisi vrednosti iz pojedinacnih zapisa
clear time u y ime_fajla n

%% Prikaz spojenih zapisa ulaza i izlaza


uneto = 0; 

while(uneto==0)
    % Da li se trazi prikaz snimljenih promena ulaza i izlaza?
    prikazi_snimljeno = input('\n Prikazati snimljeno? 0 - NE, 1 - DA!  \n');
    uneto = pravilan_unos(prikazi_snimljeno,0,1);
end

% ako je izbor takav, prikazati spojene snimljene zapise ulaza i izlaza
if (prikazi_snimljeno==1)
    
    % prikazi spojene zapise ulaza
    figure(1)
    plot(time_con,u_con, 'k')
    xlabel('$t$','Interpreter','Latex','FontSize',18)
    ylabel('$X(t)$','Interpreter','Latex','FontSize',18)

    % prikazi spojene zapise izlaza
    figure(2)
    plot(time_con,y_con, 'k')
    xlabel('$t$','Interpreter','Latex','FontSize',18)
    ylabel('$Y(t)$','Interpreter','Latex','FontSize',18)

end

%% Vrednosti wmin i wmax

% unos minimalne uchestanosti na kojoj se trazhi poznavanje uchestanosne
% karakteristike
uneto = 0; 

while(uneto==0)
    wmin = input('\n Uneti vrednost najmanje uchestanosti na kojoj se trazhi model wmin[rad/s]: \n');
    uneto = pravilan_unos(wmin,eps);
end

% ista minimalna uchestanost u [Hz]
fmin = wmin/2/pi;

% unos maksimalne uchestanosti na kojoj se trazhi poznavanje uchestanosne
% karakteristike
uneto = 0; 

while(uneto==0)
    wmax = input('\n Uneti vrednost najvece uchestanosti na kojoj se trazhi model wmax[rad/s]: \n');
    uneto = pravilan_unos(wmax,wmin);
end

% ista maksimalna uchestanost u [Hz]
fmax = wmax/2/pi;

%% Odabir broja prozora

% najmanja duzhina trajanja zapisa od svih pojedinachnih zapisa
Trecmin = min(Trec);

% najveca preporucena shirina prozora za odredjivanje spektara
Twin_max = min(Trecmin/2,TF/5);

if (Trecmin < Twin_max)
    error('\n Duzhina trajanja zapisa nije dovoljno velika. \n')
end

% najmanja preporucena shirina prozora za odredjivanje spektara
Twin_min = 2*pi*20/wmax;

%% Odabrati presechnu uchestanost filtra za sve kanale

uneto = 0;

while(uneto==0)
    wfil =  input('\n Uneti vrednost presechne uchestanosti low-pass filtra [rad/s]: \n');
    uneto = pravilan_unos(wfil,wmax,0.5*w_sampling);
end

%% Prelazak na oznake x - ulaz i y - izlaz (nista posebno)

x = u_con;
y = y_con;

%% Uklanjanje trendova iz spojenih snimaka ulaza i izlaza

y_detrend = detrend(y,'constant');
x_detrend = detrend(x,'constant');

% y_detrend = detrend(y_detrend,'linear');
% x_detrend = detrend(x_detrend,'linear');

% ako je izbor takav, prikazati spojene snimljene zapise ulaza i izlaza
% kojima su uklonjeni trendovi

uneto = 0; 

while(uneto==0)
    % Da li se trazi prikaz snimljenih promena ulaza i izlaza?
    prikazi_detrendovano = input('\n Prikazati sa uklonjenim trendovima? 0 - NE, 1 - DA!  \n');
    uneto = pravilan_unos(prikazi_detrendovano,0,1);
end

if (prikazi_detrendovano==1)

    figure(3)
    plot(time_con,y_detrend)
    xlabel('$t$','Interpreter','Latex','FontSize',18)
    ylabel('$y(t)$','Interpreter','Latex','FontSize',18)

    figure(4)
    plot(time_con,x_detrend)
    xlabel('$t$','Interpreter','Latex','FontSize',18)
    ylabel('$x(t)$','Interpreter','Latex','FontSize',18)

end

%%

s = tf('s');

%% Filtriranje snimljenih zapisa

% ako je izbor takav, prikazati spojene, isfiltrirane snimljene zapise ulaza i izlaza

uneto = 0; 

while(uneto==0)
    % Da li se trazi prikaz isfiltriranih promena ulaza i izlaza?
    prikazi_isfiltrirano = input('\n Prikazati isfiltrirane signale? 0 - NE, 1 - DA!  \n');
    uneto = pravilan_unos(prikazi_isfiltrirano,0,1);
end


% isfiltriraj zapise izlaza sa uklonjenim trendovima niskopropusnim filtrom
% drugog rada sa unetom vrednoshc1u presechne uchestanosti
[y_detrend_fil] = lsim(1/(s/(wfil)+1)/(s/(wfil)+1),y_detrend,time_con);

% isfiltriraj zapise ulaza sa uklonjenim trendovima niskopropusnim filtrom
% drugog rada sa unetom vrednoshc1u presechne uchestanosti
[x_detrend_fil] = lsim(1/(s/(wfil)+1)/(s/(wfil)+1),x_detrend,time_con);
    
if (prikazi_isfiltrirano==1)
    
%Isfiltrirati sve snimke ulaza i izlaza i prikazati uporedo

    figure(5)
    plot(time_con,y_detrend)
    hold on
    plot(time_con,y_detrend_fil,'r')
    xlabel('$t$','Interpreter','Latex','FontSize',18)
    ylabel('$y(t)$','Interpreter','Latex','FontSize',18)

    figure(6)
    plot(time_con,x_detrend)
    hold on
    plot(time_con,x_detrend_fil,'r')
    xlabel('$t$','Interpreter','Latex','FontSize',18)
    ylabel('$x(t)$','Interpreter','Latex','FontSize',18)

end

%% Odrediti spektralne vrednosti za razlichite prozore

% broj razlichitih sluchajeva 
uneto = 0; 

while(uneto==0)
    broj_razlichitih_prozora = input('\n Unesi broj razlichitih prozora koje zhelish da koristish: \n');
    uneto = pravilan_unos(broj_razlichitih_prozora,1,6);
end

% brojac do zadate vrednosti broja prozora
i = 1;

% koliko procenata od maksimalne duzhine prozora ce biti dati prozor?
perc_win = zeros(1,broj_razlichitih_prozora);

% duzhine razlichitih prozora 
Twin = zeros(1,broj_razlichitih_prozora);

% duzhine razlichitih prozora 
w_min_fit = zeros(1,broj_razlichitih_prozora);

% broj osrednjavanja u odredjenom sluchaju
nd = zeros(1,broj_razlichitih_prozora);

% izvrshi unos duzhine razlichitih prozora
% broj tachaka u prozoru
win_points = zeros(1,broj_razlichitih_prozora);

% broj tachaka koje se preklapaju
overlap_points = zeros(1,broj_razlichitih_prozora);

% procenat preklapanja prozora
uneto = 0;

while(uneto==0)
    perc_overlap = input('\n Uneti vrednost preklapanja prozora [0-1]: \n');
    uneto = pravilan_unos(perc_overlap,eps,1);
end

% unose se duzhine trajanja prozora
while (i<=broj_razlichitih_prozora)
    
    zadovoljan = 0;
    
    % dok se ne potvrdi duzhina prozora, mozhe se reci da nismo zadovoljni
    % odabirom
        
	while (zadovoljan == 0)
        
        % procenat duzhine prozora u odnosu na maksimalni dozvoljeni
        
        uneto = 0; 
        
        while(uneto==0)
            perc_win(1,i) = input('\n Unesi parametar a, koji povezuje Twin = a*Twin_max, [0-1]: \n');
            uneto = pravilan_unos(perc_win(1,i),eps,1);
        end
        
        % duzhina prozora
        Twin(1,i) = perc_win(1,i)*Twin_max;
        
        % minimalna uchestanost za fitovanje
        w_min_fit(1,i) = 2*pi/Twin(1,i);
        fprintf('\n Minimalna uchestanost za fitovanje je: %f. \n',w_min_fit(1,i));
        
        % broj osrednjavanja
        nd(1,i) = ceil(Trecmin/Twin(1,i));
        fprintf('\n Broj osrednjavanja je %f. \n',nd(1,i));
        
        if (nd(1,i) < 5)
            fprintf('\n Broj osrednjavanja je manji od 5! Probaj sa manjom vrednoscu! \n');
            break
        end
        
      	if (Twin(1,i) > Twin_max || Twin(1,i) < Twin_min )
            fprintf('\n Duzhina prozora nije odgovarajuca! \n');
            break
        end
        
        % broj tachaka u prozoru
      	win_points(1,i) = ceil((Trecmin/nd(1,i))/dT(1));

       	% broj tachaka koje se preklapaju
      	overlap_points(1,i) = ceil(perc_overlap*(Trecmin/nd(1,i))/dT(1));
        
      	% potvrditi broj osrednjvanja
      	zadovoljan = input('\n Unesi: 1 - zadovoljan, 0 - nezadovoljan!  \n');
        
       	% proveriti da li su uneti 1 ili 0
      	if (zadovoljan == 1)
            
            i = i + 1;
            
        elseif (zadovoljan ~= 0)
            
            fprintf('\n Unos nije sadrzhao 0 ili 1. \n');
            break;
            
       	end
            
	end
        
end

clear i 


%% inicijalizacija matrica za skladishtenje srachunatih vrednosti

% matrica autospektara ulaza za razlichite prozore
Gxx = zeros(2^(nextpow2(length(x_detrend)))/2+1,broj_razlichitih_prozora);

% matrica autospektara izlaza za razlichite prozore
Gyy = zeros(2^(nextpow2(length(y_detrend)))/2+1,broj_razlichitih_prozora);

% matrica medjusobnih spekatara ulaza i izlaza za razlichite prozore
Gxy = zeros(2^(nextpow2(length(x_detrend)))/2+1,broj_razlichitih_prozora);
Gyx = zeros(2^(nextpow2(length(y_detrend)))/2+1,broj_razlichitih_prozora);

% matrica funkcija koherencije za razlichite prozore
gamma_xy = zeros(2^(nextpow2(length(x_detrend)))/2+1,broj_razlichitih_prozora);

% matrica funkcija nasumichne greshke
epsr = zeros(2^(nextpow2(length(x_detrend)))/2+1,broj_razlichitih_prozora);

% matrica vrednosti uchestanosnih karakteristika
H = zeros(2^(nextpow2(length(x_detrend)))/2+1,broj_razlichitih_prozora);

% matrica vrednosti log-amp uchestanosnih karakteristika
L_Amp = zeros(2^(nextpow2(length(x_detrend)))/2+1,broj_razlichitih_prozora);

% matrica vrednosti faznih uchestanosnih karakteristika
Phi = zeros(2^(nextpow2(length(x_detrend)))/2+1,broj_razlichitih_prozora);

%% Srachunati trazhene velichine

% brojach vezan za nijansu boje prikaza
k = 0;

for i = 1:broj_razlichitih_prozora
    
% Izvrshiti estimaciju auto- i medjusobnih spektara

    % auto-spektar ulaza 
    [Gxx(:,i),f] = cpsd(x_detrend_fil,x_detrend_fil,...
                        hanning(win_points(1,i)),overlap_points(1,i),...
                        2^(nextpow2(length(x_detrend))),f_sampling);                

    figure(7)
    semilogx(f*2*pi,10*log10(Gxx(:,i)),'Color',[0/255 k/255 255/255])
    hold on
    
    if (i==broj_razlichitih_prozora)

        i1 = find(2*pi*f<=0.8*wmin,1,'last');
        i2 = find(2*pi*f>=1.2*wmax,1,'first');

        Ylow = min(min( 10*log10(Gxx(i1:i2,:)) ));
        Yhigh = max(max( 10*log10(Gxx(i1:i2,:)) ));

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

        xlabel('$\omega [\mathrm{rad/s}]$','Interpreter',...
        'Latex','FontSize',18);
        ylabel('$G_{xx}(\omega)[\mathrm{dB}]$','Interpreter',...
        'Latex','FontSize',18);
        set(gcf,'Position',[10 10 800 450])
        
    end

    % auto-spektar izlaza
    [Gyy(:,i),f] = cpsd(y_detrend_fil,y_detrend_fil,...
                        hanning(win_points(1,i)),overlap_points(1,i),...
                        2^(nextpow2(length(y_detrend))),f_sampling);

    figure(8)
    semilogx(f*2*pi,10*log10(Gyy(:,i)),'Color',[0/255 k/255 255/255])
    hold on
    
    if (i==broj_razlichitih_prozora)

        i1 = find(2*pi*f<=0.8*wmin,1,'last');
        i2 = find(2*pi*f>=1.2*wmax,1,'first');

        Ylow = min(min( 10*log10(Gyy(i1:i2,:)) ));
        Yhigh = max(max( 10*log10(Gyy(i1:i2,:)) ));

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

        xlabel('$\omega [\mathrm{rad/s}]$','Interpreter',...
        'Latex','FontSize',18);
        ylabel('$G_{yy}(\omega)[\mathrm{dB}]$','Interpreter',...
        'Latex','FontSize',18);
        set(gcf,'Position',[10 10 800 450])

    end

    % Medjusobni spektri ulaza i izlaza
    [Gxy(:,i),~] = cpsd(y_detrend_fil,x_detrend_fil,...
                        hanning(win_points(1,i)),overlap_points(1,i),...
                        2^(nextpow2(length(x_detrend))),f_sampling);
    
    figure(18)
    subplot(2,1,1)
    semilogx(f*2*pi,10*log10(abs(Gxy(:,i))),'Color',[0/255 k/255 255/255])
    hold on
    if (i==broj_razlichitih_prozora)

        i1 = find(2*pi*f<=0.8*wmin,1,'last');
        i2 = find(2*pi*f>=1.2*wmax,1,'first');

        Ylow = min(min( 10*log10(abs(Gxy(i1:i2,:))) ));
        Yhigh = max(max( 10*log10(abs(Gxy(i1:i2,:))) ));

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

        xlabel('$\omega [\mathrm{rad/s}]$','Interpreter',...
        'Latex','FontSize',16);
        ylabel('$|G_{xy}|(\omega)[\mathrm{dB}]$','Interpreter',...
        'Latex','FontSize',16);
        %set(gcf,'Position',[10 10 800 450])
    end
    subplot(2,1,2)
    semilogx(f*2*pi,angle(Gxy(:,i)),'b','Color',[0/255 k/255 255/255])
    hold on
    if (i==broj_razlichitih_prozora)
        
        i1 = find(2*pi*f<=0.8*wmin,1,'last');
        i2 = find(2*pi*f>=1.2*wmax,1,'first');
        
        Ylow = min(min( angle(Gxy(i1:i2,:)) ));
        Yhigh = max(max( angle(Gxy(i1:i2,:)) ));
        
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
        
        xlabel('$\omega [\mathrm{rad/s}]$','Interpreter',...
            'Latex','FontSize',16);
        ylabel('$ \angle (G_{xy})(\omega) $','Interpreter',...
            'Latex','FontSize',16);
        set(gcf,'Position',[10 -200 1000 600])
        
    end
    
    
    [Gyx(:,i),~] = cpsd(x_detrend_fil,y_detrend_fil,...
                        hanning(win_points(1,i)),overlap_points(1,i),...
                        2^(nextpow2(length(x_detrend))),f_sampling);
    
    figure(19)
    subplot(2,1,1)
    semilogx(f*2*pi,10*log10(abs(Gyx(:,i))),'Color',[0/255 k/255 255/255])
    hold on
    if (i==broj_razlichitih_prozora)

        i1 = find(2*pi*f<=0.8*wmin,1,'last');
        i2 = find(2*pi*f>=1.2*wmax,1,'first');

        Ylow = min(min( 10*log10(abs(Gyx(i1:i2,:))) ));
        Yhigh = max(max( 10*log10(abs(Gyx(i1:i2,:))) ));

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

        xlabel('$\omega [\mathrm{rad/s}]$','Interpreter',...
        'Latex','FontSize',16);
        ylabel('$|G_{yx}|(\omega)[\mathrm{dB}]$','Interpreter',...
        'Latex','FontSize',16);
        %set(gcf,'Position',[10 10 800 450])
    end
    subplot(2,1,2)
    semilogx(f*2*pi,angle(Gyx(:,i)),'b','Color',[0/255 k/255 255/255])
    hold on
    if (i==broj_razlichitih_prozora)
        
        i1 = find(2*pi*f<=0.8*wmin,1,'last');
        i2 = find(2*pi*f>=1.2*wmax,1,'first');
        
        Ylow = min(min( angle(Gyx(i1:i2,:)) ));
        Yhigh = max(max( angle(Gyx(i1:i2,:)) ));
        
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
        
        xlabel('$\omega [\mathrm{rad/s}]$','Interpreter',...
            'Latex','FontSize',16);
        ylabel('$ \angle (G_{yx})(\omega) $','Interpreter',...
            'Latex','FontSize',16);
        set(gcf,'Position',[10 -200 1000 600])
        
    end
        
    

    % Procene uchestanosnih karakteristika
    H1 = Gxy(:,i)./Gxx(:,i);
    H2 = Gyy(:,i)./Gyx(:,i);

    % optimalno resenje za slucaj kada nema shuma na ulazu
    H(:,i) = H1;
    L_Amp(:,i) = 20*log10(abs(H(:,i)));
    Phi(:,i) = unwrap(angle(H(:,i)))*180/pi;

% Funkcija koherencije

    figure(9)
    [gamma_xy(:,i),~] = mscohere(x_detrend_fil,y_detrend_fil,...
                                 hanning(win_points(1,i)),overlap_points(1,i),...
                                 2^(nextpow2(length(x_detrend_fil))),f_sampling);
    semilogx(f*2*pi,gamma_xy(:,i),'Color',[255/255 k/255 0])
    hold on
    
    if (i==broj_razlichitih_prozora)

        i1 = find(2*pi*f<=0.8*wmin,1,'last');
        i2 = find(2*pi*f>=1.2*wmax,1,'first');

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
        hold on

        xlabel('$\omega [\mathrm{rad/s}]$','Interpreter',...
        'Latex','FontSize',18);
        ylabel('$\gamma_{xy}^2(\omega)[\mathrm{dB}]$','Interpreter',...
        'Latex','FontSize',18);
        set(gcf,'Position',[10 10 800 450])
        legend({'$a = 0.25$', '$a = 0.3$', '$ a = 0.35$'}, 'Interpreter','Latex', 'Location', 'southwest')
        
    end

    % Funkcija greshke

    % funkcija nasumichne greshke
    epsr(:,i) = sqrt(0.5)*sqrt(1-gamma_xy(:,i))./abs(gamma_xy(:,i))./sqrt(2*nd(1,i));

    % Prikaz procenjene uchestanosne karakteristike

    figure(10)
    subplot(3,1,1)
    semilogx(f*2*pi,L_Amp(:,i),'Color',[0/255 k/255 255/255])
    hold on
    
    if (i==broj_razlichitih_prozora)

        i1 = find(2*pi*f<=0.8*wmin,1,'last');
        i2 = find(2*pi*f>=1.2*wmax,1,'first');

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

        xlabel('$\omega [\mathrm{rad/s}]$','Interpreter',...
               'Latex','FontSize',14);
        ylabel('$ L(\omega)[\mathrm{dB}] $','Interpreter',...
               'Latex','FontSize',14);
        legend({'$a = 0.25$', '$a = 0.3$', '$ a = 0.35$'}, 'Interpreter','Latex')
    end

    subplot(3,1,2)
    semilogx(f*2*pi,Phi(:,i),'b','Color',[0/255 k/255 255/255])
    hold on
    
    if (i==broj_razlichitih_prozora)

        i1 = find(2*pi*f<=0.8*wmin,1,'last');
        i2 = find(2*pi*f>=1.2*wmax,1,'first');

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

        xlabel('$\omega [\mathrm{rad/s}]$','Interpreter',...
               'Latex','FontSize',14);
        ylabel('$ \varphi(\omega) [^{\circ}]$','Interpreter',...
               'Latex','FontSize',14);
           legend({'$a = 0.25$', '$a = 0.3$', '$ a = 0.35$'}, 'Interpreter','Latex', 'Location', 'southwest')

    end

    subplot(3,1,3)
    semilogx(f*2*pi,gamma_xy(:,i),'Color',[255/255 k/255 0])
    hold on
    
    if (i==broj_razlichitih_prozora)

        i1 = find(2*pi*f<=0.8*wmin,1,'last');
        i2 = find(2*pi*f>=1.2*wmax,1,'first');

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

        xlabel('$\omega [\mathrm{rad/s}]$','Interpreter',...
               'Latex','FontSize',14);
        ylabel('$ \gamma_xy^2(\omega) $','Interpreter',...
               'Latex','FontSize',14);
        semilogx([0.8*wmin 1.2*wmax],[0.6 0.6],'k--')
        set(gcf,'Position',[10 -200 800 1200])
        legend({'$a = 0.25$', '$a = 0.3$', '$ a = 0.35$'}, 'Interpreter','Latex', 'Location', 'southwest')
        
    end
    
    k = k+ 255/(broj_razlichitih_prozora);
    
end

%% Vektor uchestanosti u [rad/s]

w = 2*f*pi;

%% Snimanje rezultata rada

dT = dT(1);

ime_fajla = input('\n Unesi ime eksperimenta: \n');

save([ime_fajla '_frespid' '.mat'],'broj_razlichitih_prozora','dT','wmin','wmax',...
                        'w','w_min_fit','Twin','nd','overlap_points',...
                        'win_points','x_detrend','y_detrend','Gxx','Gyy',...
                        'Gxy','Gyx','gamma_xy','epsr','H','L_Amp','Phi')