hold off
Vmin  = 0.0001;                 % depth of potential
alpha = 3.5000;
theta = [Vmin alpha];

a1 = 0.004;
a2 = 0.00005;
k1 = 10;
k2 = 2;
l1 = 0.40;
l2 = 1.6;
greta = [k1 l1 a1 k2 l2 a2];

r = linspace(0,5,1001);
y = U_pot(r,theta);
z = U_pot(r,greta);

plot(r,y,'b')
hold on
plot(r,z,'r')
grid on
axis([0 5 -0.001 0.005])