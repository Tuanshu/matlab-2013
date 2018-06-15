clear all

cd('D:\Users\TuanShu\');

Thickness_substrate=500;    %micron
Thickness_film=4.5;         %micron

n_film=importdata('n_5micron_NewAlgorithm.txt');
Wavelength_micron_Considered=importdata('Wavelength_micron_Considered_5micron_NewAlgorithm.txt');

Wavelength_center=0.560;    %micron

c=3E8;

NA_Obj=0.25;

DOF=Wavelength_center/NA_Obj^2;

Depth=(-100+0.001):0.001:100;

PSF_DOF=gaussmf(Depth,[DOF 0]);

%% Substrate RI

% Assumed n1 = BK7
C1 = 1.03961212; 
C2 = 0.00600069867; 
C3 = 0.231792344; 
C4 = 0.0200179144; 
C5 = 1.01046945; 
C6 = 103.560653;


n_bk7=(C1*(Wavelength_micron_Considered.^2)./((Wavelength_micron_Considered.^2)-C2)+C3*(Wavelength_micron_Considered.^2)./((Wavelength_micron_Considered.^2)-C4)+C5*(Wavelength_micron_Considered.^2)./((Wavelength_micron_Considered.^2)-C6)+1).^0.5;

n_bk7=abs(n_bk7);

n_bk7(isnan(n_bk7))=0;
n_substrate=n_bk7;%/1.516781257666726*1.520;

%% Defocus

Defocus_substrate=Thickness_substrate.*(n_substrate-1./n_substrate);
Defocus_substrate=Defocus_substrate-Defocus_substrate(find(Wavelength_micron_Considered<Wavelength_center,1,'first'));

for p=1:length(Defocus_substrate)
    Ratio_substrate(p)=PSF_DOF(find(Depth>Defocus_substrate(p),1,'first'));
end


Defocus_film=Thickness_film.*(n_film-1./n_film);
Defocus_film=Defocus_film-Defocus_film(find(Wavelength_micron_Considered<Wavelength_center,1,'first'));

for p=1:length(Defocus_film)
    Ratio_film(p)=PSF_DOF(find(Depth>Defocus_film(p),1,'first'));
end


plot(Wavelength_micron_Considered,Defocus_substrate);
xlabel('Wavelength (micron)');
ylabel('Defocus (micron)');


plot(Wavelength_micron_Considered,Defocus_film);
xlabel('Wavelength (micron)');
ylabel('Defocus (micron)');


plot(Wavelength_micron_Considered,n_film);
xlabel('Wavelength (micron)');
ylabel('n');


plot(Wavelength_micron_Considered,Defocus_substrate,Wavelength_micron_Considered,Defocus_film);

plot(Wavelength_micron_Considered,10*log10(Ratio_substrate));
xlabel('Wavelength (micron)');
ylabel('Loss (dB)');

plot(Wavelength_micron_Considered,10*log10(Ratio_film));
xlabel('Wavelength (micron)');
ylabel('Loss (dB)');
