clear all

Spectra_Array=[10 15 20 25 30 35 40];
%Spectra_Array=[5];

Path='D:\Users\TuanShu\130225_CeLYSO Spectra';

Max_Wavelength=550;
Min_Wavelength=360;

Extra_Amp_Coef_of_3rd_Term=1;

ToleranceTG=1E-6;
NofIerationTG=1000;

cd(sprintf('%s\\',Path));

for p=1:length(Spectra_Array)
    Data=dlmread(sprintf('%dmW.txt',Spectra_Array(p)));
    if p==1
        Wavelength=Data(:,1);   %nm
        Spectra(1:length(Wavelength),1:length(Spectra_Array))=0;
    end
    Spectra(:,p)=Data(:,3);
end
    
df=abs(3E8/(Wavelength(length(Wavelength))*1E-9)-3E8/(Wavelength(length(Wavelength)-1)*1E-9));

%plot(Wavelength,Spectra);

Frequency=3E8./(Wavelength*1E-9);

Frequency_New=(3E8/(Max_Wavelength*1E-9)):df:(3E8/(Min_Wavelength*1E-9));
Frequency_New=Frequency_New';
Spectra_New=interp1(Frequency,Spectra,Frequency_New);
Spectra_New(isnan(Spectra_New))=0;
Frequency_New=Frequency_New/1E14;

Starting_Frequency_1= 7.1481;  
Starting_Frequency_2=7.6876;
Starting_Frequency_3=6.3881;
Starting_dFrequency_1=0.3619;
Starting_dFrequency_2=0.2091;
Starting_dFrequency_3=0.3223;
Starting_Amp_1=0.8898;
Starting_Amp_2=0.6538;
Starting_Amp_3=0.113;    %first assume only the amp changes with power

%FitDualGaussian= fittype( @(Parameter_Frequency_1,Parameter_Frequency_2,Parameter_dFrequency_1,Parameter_dFrequency_2,Parameter_Amp_1,Parameter_Amp_2,f) Parameter_Amp_1*gaussmf(f,[Parameter_dFrequency_1 Parameter_Frequency_1])+Parameter_Amp_2*gaussmf(f,[Parameter_dFrequency_2 Parameter_Frequency_2]), 'independent', {'f'});%,'dependent', 'z');    %(x-A)^2+(y-B)^2+(z-C)^2=R^2
FitTriGaussian=fittype( @(Parameter_Frequency_1,Parameter_Frequency_2,Parameter_Frequency_3,Parameter_dFrequency_1,Parameter_dFrequency_2,Parameter_dFrequency_3,Parameter_Amp_1,Parameter_Amp_2,Parameter_Amp_3,f) Parameter_Amp_1*gaussmf(f,[Parameter_dFrequency_1 Parameter_Frequency_1])+Parameter_Amp_2*gaussmf(f,[Parameter_dFrequency_2 Parameter_Frequency_2])+Extra_Amp_Coef_of_3rd_Term*Parameter_Amp_3*gaussmf(f,[Parameter_dFrequency_3 Parameter_Frequency_3]), 'independent', {'f'});%,'dependent', 'z');    %(x-A)^2+(y-B)^2+(z-C)^2=R^2

%FitOptionsForTG=fitoptions('Method','NonlinearLeastSquares','StartPoint',[Starting_Frequency_1,Starting_Frequency_2,Starting_Frequency_3,Starting_dFrequency_1,Starting_dFrequency_2,Starting_dFrequency_3,Starting_Amp_1,Starting_Amp_2,Starting_Amp_3],'TolFun',ToleranceTG,'MaxIter',NofIerationTG);
FitOptionsForTG=fitoptions('Method','NonlinearLeastSquares','StartPoint',[Starting_Frequency_1,Starting_Frequency_2,Starting_Frequency_3,Starting_dFrequency_1,Starting_dFrequency_2,Starting_dFrequency_3,Starting_Amp_1,Starting_Amp_2,Starting_Amp_3],'TolFun',ToleranceTG,'MaxIter',NofIerationTG);

StartCurve=Starting_Amp_1*gaussmf(Frequency_New,[Starting_dFrequency_1 Starting_Frequency_1])+Starting_Amp_2*gaussmf(Frequency_New,[Starting_dFrequency_2 Starting_Frequency_2])+Extra_Amp_Coef_of_3rd_Term*Starting_Amp_3*gaussmf(Frequency_New,[Starting_dFrequency_3 Starting_Frequency_3]);
CheckCurve(1:length(Frequency_New),1:length(Spectra_Array))=0;
Parameter_Array_Frequency_1(1:length(Spectra_Array))=0;
Parameter_Array_Frequency_2(1:length(Spectra_Array))=0;
Parameter_Array_Frequency_3(1:length(Spectra_Array))=0;
Parameter_Array_dFrequency_1(1:length(Spectra_Array))=0;
Parameter_Array_dFrequency_2(1:length(Spectra_Array))=0;
Parameter_Array_dFrequency_3(1:length(Spectra_Array))=0;
Parameter_Array_Amp_1(1:length(Spectra_Array))=0;
Parameter_Array_Amp_2(1:length(Spectra_Array))=0;
Parameter_Array_Amp_3(1:length(Spectra_Array))=0;
ErrorCalc(1:length(Spectra_Array))=0;
Parameter_Array_Frequency_1=Parameter_Array_Frequency_1';
Parameter_Array_Frequency_2=Parameter_Array_Frequency_2';
Parameter_Array_Frequency_3=Parameter_Array_Frequency_3';
Parameter_Array_dFrequency_1=Parameter_Array_dFrequency_1';
Parameter_Array_dFrequency_2=Parameter_Array_dFrequency_2';
Parameter_Array_dFrequency_3=Parameter_Array_dFrequency_3';
Parameter_Array_Amp_1=Parameter_Array_Amp_1';
Parameter_Array_Amp_2=Parameter_Array_Amp_2';
Parameter_Array_Amp_3=Parameter_Array_Amp_3';
ErrorCalc=ErrorCalc';
for q=1:length(Spectra_Array)
    [FitResult Error]=fit([Frequency_New],Spectra_New(:,q),FitTriGaussian,FitOptionsForTG);
    CheckCurve(:,q)=FitResult.Parameter_Amp_1*gaussmf(Frequency_New,[FitResult.Parameter_dFrequency_1 FitResult.Parameter_Frequency_1])+FitResult.Parameter_Amp_2*gaussmf(Frequency_New,[FitResult.Parameter_dFrequency_2 FitResult.Parameter_Frequency_2])+Extra_Amp_Coef_of_3rd_Term*FitResult.Parameter_Amp_3*gaussmf(Frequency_New,[FitResult.Parameter_dFrequency_3 FitResult.Parameter_Frequency_3]);
    Parameter_Array_Frequency_1(q)=FitResult.Parameter_Frequency_1;
    Parameter_Array_Frequency_2(q)=FitResult.Parameter_Frequency_2;
    Parameter_Array_Frequency_3(q)=FitResult.Parameter_Frequency_3;
    Parameter_Array_dFrequency_1(q)=FitResult.Parameter_dFrequency_1;
    Parameter_Array_dFrequency_2(q)=FitResult.Parameter_dFrequency_2;
    Parameter_Array_dFrequency_3(q)=FitResult.Parameter_dFrequency_3;
    Parameter_Array_Amp_1(q)=FitResult.Parameter_Amp_1;
    Parameter_Array_Amp_2(q)=FitResult.Parameter_Amp_2;
    Parameter_Array_Amp_3(q)=FitResult.Parameter_Amp_3;
    ErrorCalc(q)=sum((CheckCurve(:,q)-Spectra_New(:,q)).^2);
    disp(q);
end
Array_Frequency=[Parameter_Array_Frequency_1 Parameter_Array_Frequency_2 Parameter_Array_Frequency_3];
Array_dFrequency=[Parameter_Array_dFrequency_1 Parameter_Array_dFrequency_2 Parameter_Array_dFrequency_3];
Array_Amp=[Parameter_Array_Amp_1 Parameter_Array_Amp_2 Parameter_Array_Amp_3];

Array_Wavelength=(3E8./(Array_Frequency*1E14))*1E9;  %nm           df/f=dl/l
Array_dWavelength=Array_Wavelength.*(Array_dFrequency./Array_Frequency);
%%
Show_Number=1;
plot(Frequency_New,CheckCurve(:,Show_Number),Frequency_New,Spectra_New(:,Show_Number),Frequency_New,StartCurve);
%plot(Frequency_New,CheckCurve(:,Show_Number));
Parameter_Array_Frequency_1(Show_Number)
Parameter_Array_Frequency_2(Show_Number)
Parameter_Array_Frequency_3(Show_Number)
Parameter_Array_dFrequency_1(Show_Number)
Parameter_Array_dFrequency_2(Show_Number)
Parameter_Array_dFrequency_3(Show_Number)
Parameter_Array_Amp_1(Show_Number)
Parameter_Array_Amp_2(Show_Number)
Parameter_Array_Amp_3(Show_Number)
ErrorCalc(Show_Number)

plot(Spectra_Array,ErrorCalc);
xlabel('Pump Power (mW)');
ylabel('Error (std)');

plot(Spectra_Array,Array_Wavelength);
xlabel('Pump Power (mW)');
ylabel('Center Wavelength of Gaussian (nm)');

plot(Spectra_Array,Array_dWavelength);
xlabel('Pump Power (mW)');
ylabel('FWHM of Gaussian (nm)');

plot(Spectra_Array,Array_dWavelength);
xlabel('Pump Power (mW)');
ylabel('Amplitude of Gaussian (nm)');


dlmwrite('130225_Error_3 Gauss_All Variable.txt',ErrorCalc);
dlmwrite('130225_Center Frequency_3 Gauss_All Variable.txt',Array_Frequency);
dlmwrite('130225_dFrequency_3 Gauss_All Variable.txt',Array_Frequency);
dlmwrite('130225_Center Wavelength_3 Gauss_All Variable.txt',Array_Wavelength);
dlmwrite('130225_dWavelength_3 Gauss_All Variable.txt',Array_dWavelength);
dlmwrite('130225_dWavelength_3 Gauss_All Variable.txt',Array_dWavelength);
dlmwrite('130225_Amplitude_3 Gauss_All Variable.txt',Array_Amp);

%Error