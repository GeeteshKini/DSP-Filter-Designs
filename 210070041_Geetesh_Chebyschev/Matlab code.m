r = [-0.1125-0.79445j -0.1125+0.79445j -0.04015-0.99067j -0.04015+0.99067j -0.16255+0.44089j -0.16255-0.44089j -0.18041]
p = poly(r)

%Chebyshev LPF parameters
D1 = 1/(0.85*0.85)-1;       %since delta is 0.15
epsilon = sqrt(D1);         %epsilon was set to this value to satisfy required inequality
N = 7;

%Parameters for Bandpass Transformation
W0 = 0.61945;
B = 0.573;

%Evaluating Frequency Response of Final Filter
syms s z;
num=[0.0252];
den= @(s) 1*s.^7 + 0.81*s.^6 + 2.0785*s.^5 + 1.2284*s.^4 + 1.2427*s.^3 + 0.4613*s.^2 + 0.1878*s + 0.0252
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);    %analog lpf transfer function
analog_bpf(s) = analog_lpf((s*s +W0*W0)/(B*s));     %bandpass transformation
discrete_bpf(z) = analog_bpf((z-1)/(z+1));          %bilinear transformation

%coeffs of analog BPF
[ns, ds] = numden(analog_bpf(s));                   %numerical simplification
ns = sym2poly(expand(ns));                          
ds = sym2poly(expand(ds));                          %collect coeffs into matrix form
k = ds(1);    
ds = ds/k;
ns = ns/k;
ds
ns

%coeffs of discrete BPF
[nz, dz] = numden(discrete_bpf(z));                 %numerical simplification
nz = sym2poly(expand(nz));                          
dz = sym2poly(expand(dz));                          %collect coeffs into matrix form
k = dz(1);                                          %normalisation factor
dz = dz/k;
nz = nz/k;
dz
nz

fvtool(nz,dz)
%magnitude plot (not in log scale) 
[H,f] = freqz(nz,dz,1024*1024, 600e3);
plot(f,abs(H))
plot(f,angle(H))
grid