clear all

Path='D:\Users\TuanShu\FWHM.txt';

Data=importdata(Path);      % Data_2: the glass data 1

Data=Data';


dlmwrite(Path,Data,'delimiter','\t','newline','pc');

