%Butterworth Analog LPF parameters
Wc = 1.05;              %cut-off frequency
N = 12;                  %order 

%poles of Butterworth polynomial of degree 12 in the open CLHP 
p1 = Wc*cos(pi/2 + pi/24) + i*Wc*sin(pi/2 + pi/24);
p2 = Wc*cos(pi/2 + pi/24) - i*Wc*sin(pi/2 + pi/24);
p3 = Wc*cos(pi/2 + pi/24+pi/12) + i*Wc*sin(pi/2 + pi/24+pi/12);
p4 = Wc*cos(pi/2 + pi/24+pi/12) - i*Wc*sin(pi/2 + pi/24+pi/12);
p5 = Wc*cos(pi/2 + pi/24+2*pi/12) + i*Wc*sin(pi/2 + pi/24+2*pi/12);
p6 = Wc*cos(pi/2 + pi/24+2*pi/12) - i*Wc*sin(pi/2 + pi/24+2*pi/12);
p7 = Wc*cos(pi/2 + pi/24+3*pi/12) + i*Wc*sin(pi/2 + pi/24+3*pi/12);
p8 = Wc*cos(pi/2 + pi/24+3*pi/12) - i*Wc*sin(pi/2 + pi/24+3*pi/12);
p9 = Wc*cos(pi/2 + pi/24+4*pi/12) + i*Wc*sin(pi/2 + pi/24+4*pi/12);
p10 = Wc*cos(pi/2 + pi/24+4*pi/12) - i*Wc*sin(pi/2 + pi/24+4*pi/12);
p11 = Wc*cos(pi/2 + pi/24+5*pi/12) + i*Wc*sin(pi/2 + pi/24+5*pi/12);
p12 = Wc*cos(pi/2 + pi/24+5*pi/12) - i*Wc*sin(pi/2 + pi/24+5*pi/12);

%Band Edge speifications
fp1 = 65;
fs1 = 70;
fs2 = 110;
fp2 = 115;

%Transformed Band Edge specs using Bilinear Transformation
f_samp = 425;         
wp1 = tan(fp1/f_samp*pi);
ws1 = tan(fs1/f_samp*pi); 
ws2 = tan(fs2/f_samp*pi);
wp2 = tan(fp2/f_samp*pi);

%Parameters for Bandstop Transformation
W0 = sqrt(wp1*wp2);
B = wp2-wp1;

[num,den] = zp2tf([],[p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12],Wc^N);   %TF with poles p1-p12 and numerator Wc^N and no zeroes
                                                        %numerator chosen to make the DC Gain = 1

%Evaluating Frequency Response of Final Filter
syms s z;
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);        %analog LPF Transfer Function
analog_bsf(s) = analog_lpf((B*s)/(s*s + W0*W0));        %bandstop transformation
discrete_bsf(z) = analog_bsf((z-1)/(z+1));              %bilinear transformation

%coeffs of analog bsf
[ns, ds] = numden(analog_bsf(s));                   %numerical simplification to collect coeffs
ns = sym2poly(expand(ns));                          
ds = sym2poly(expand(ds));                          %collect coeffs into matrix form
k = ds(1);    
ds = ds/k;
ns = ns/k;

%coeffs of discrete bsf
[nz, dz] = numden(discrete_bsf(z));                     %numerical simplification to collect coeffs                    
nz = sym2poly(expand(nz));
dz = sym2poly(expand(dz));                              %coeffs to matrix form
k = dz(1);                                              %normalisation factor
dz = dz/k;
nz = nz/k;
fvtool(nz,dz)                                           %frequency response

%magnitude plot (not in log scale) 
[H,f] = freqz(nz,dz,1024*1024, 425e3);
plot(f,abs(H))
grid