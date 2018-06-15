clear all

% Light source
NA_s=0.6;
w_s=40;     %micron

% lenses
f_IOL=18;
f_OBJ=4.6;  %mm
f_CL1=50;
f_CL2=50;
f_SL=18;
f_EP=21;
f_CL3=50;
f_CL4=50;
f_L1=50;
f_L2=100;
f_L3=90;

% Features

feature_FOL=2*NA_s*(f_OBJ/f_CL1)*(f_SL/f_EP)*f_IOL/1.33;          %mm
feature_M=(f_EP/f_IOL)*(f_CL3/f_SL)*(f_L1/f_CL4)*(f_L3/f_L2);
feature_b_camera=feature_FOL*feature_M;                                         %mm
feature_w_camera=w_s*(f_L1/f_OBJ)*(f_L3/f_L2);                  %micron
feature_w_retina=w_s*(f_SL/f_OBJ)*(f_IOL/f_EP);                 %micron
feature_b_grating_h=2*NA_s*f_OBJ*(f_L2/f_L1);                     %mm