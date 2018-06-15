clear all;

pitch=1000/300;        %micron, grating pitch
center_wavelength=0.76;     %micron, center wavelength
incidence_angle=8;      %degree
s_wavelength=0.6;           %shortest wavelength cared
l_wavelength=1;           %longest wavelength cared

Camera_Size=11;             %mm


s_angle=asin(s_wavelength./pitch-sin(incidence_angle/180*pi))*180/pi;   %degree
l_angle=asin(l_wavelength./pitch-sin(incidence_angle/180*pi))*180/pi;   %degree
full_angle=l_angle-s_angle;
blazed_angle=asin(center_wavelength/2/pitch)*180/pi;   %degree

required_focus=Camera_Size/2/tan(full_angle/180*pi/2);