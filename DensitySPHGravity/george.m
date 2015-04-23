function [x] = george(fn)
l = load(fn);
r = l(:,2);
a = l(:,5);
rho = l(:,3);

figure
subplot(2, 3, 1)
hold on
loglog(r,rho,'.')
ra = [1e-2 : 1e-2 : 1000];
rhoa = 3/(4*3.14) * 1./(1+ra.^2).^(5/2);
loglog(ra,rhoa,'r.');
xlabel("r")
ylabel('\rho')

subplot(2, 3, 2)
hold on
aa = ra./(1 + ra.^2).^(3/2);
loglog(ra,aa,'r.')
loglog(r,a,'.')
xlabel("r")
ylabel('a')

subplot(2, 3, 3)
hold on
phi = l(:,6);
phia = 1./(1. + ra.^2).^(1/2);
loglog(ra,phia,'r.');
loglog(r,-phi,'.');
 xlabel("r")
 ylabel('\phi')

subplot(2, 3, 4)
 nn = l(:,4);
 semilogx(r,nn,'.')
 xlabel("r")
 ylabel('nn')

subplot(2, 3, 5)
 hist(log10(nn), 20)
 xlabel('nn')
 ylabel('#')

subplot(2, 3, 6)
 h = l(:,7);
 loglog(h, r, '.')
 xlabel("r")
 ylabel('h')

end
