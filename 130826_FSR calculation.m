clear all

C=3E7;
Center_wavelength=1.064;    %micron
Cavity_length=15;            %mm
FSR=C/Cavity_length;
Center_Frequency=C/(Center_wavelength*1E-6);
delta_wavelength=Center_wavelength*FSR/Center_Frequency;
