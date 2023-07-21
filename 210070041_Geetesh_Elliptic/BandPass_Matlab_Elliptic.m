%Band Edge speifications
fp1 = 72;
fs1 = 67;
fs2 = 152;
fp2 = 147;

%Transformed Band Edge specs using Bilinear Transformation
f_samp = 600;         
wp1 = tan(fp1/f_samp*pi);
ws1 = tan(fs1/f_samp*pi); 
ws2 = tan(fs2/f_samp*pi);
wp2 = tan(fp2/f_samp*pi);

%Parameters for Bandstop Transformation
w0 = sqrt(wp1*wp2);
B = wp2-wp1;
del = 0.15;
omega_p = 1;
omega_s = 1.128;

epsilon_p = 0.6197;
epsilon_s = 6.5912;

a1= epsilon_p/epsilon_s;
a = omega_p/omega_s;

syms phi k F(phi,k) K(k) Kprime(k) 
F(phi,k) = 1./sqrt(1-(k.*sin(phi)).^2);

K(k) = int(F(phi,k),phi,0,pi/2);
Kprime(k) = subs(K(k),k,sqrt(1-k^2));

N = eval((feval(K,a)*feval(Kprime,a1))/(feval(K,a1)*feval(Kprime,a)));
N = ceil(N);

a = double(ellipdeg(N,a1));

W_snew = 1/a;

l = floor(N/2);
r = mod(N,2);

i = 1:1:l;
u = double(((2*i)-1)/N);
zeta = double(cde(u,a));
h_zeroes = double((1j)./(a.*zeta));
h_zeroes = [h_zeroes,conj(h_zeroes(1:l))];

syms x y darcsn(x,y) arcsn(x,y)
darcsn(x,y) = 1./(sqrt(1-x.^2)*sqrt(1-(x*y).^2))
arcsn(x,y) = int(darcsn(x,y),x,0,x);

v = double(-1j*asne(1i/epsilon_p,a1)/N);

h_poles = 1j*double(cde((u-1j*v),a));

if(r==1)
    p0 = 1j*double(cde((1-1j*v),a));
    h_poles = [h_poles,conj(h_poles(1:l)),p0];
else
    h_poles = [h_poles,conj(h_poles(1:l))];
end

h_poles
num = 0.15*poly(h_zeroes)
den = real(poly(h_poles))

syms s

analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);
analog_bsf(s) = analog_lpf((s*s+w0*w0)/(B*s));
[ns,ds] = numden(analog_bsf(s));
ns = sym2poly(expand(ns));
ds = sym2poly(expand(ds));
k = ds(1);
ns = ns/k
ds = ds/k

syms z
discrete_bsf(z) = analog_bsf((z-1)/(z+1));
[nz,dz] = numden(discrete_bsf(z));
nz = sym2poly(expand(nz));
dz = sym2poly(expand(dz));
k = dz(1);
nz = nz/k
dz = dz/k

fvtool(nz,dz)

[H,f] = freqz(nz,dz,1024*1024,600e3);
plot(f,abs(H))
 