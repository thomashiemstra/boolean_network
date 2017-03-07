test = figure( 'Name', 'DoS N=256');
semilogx(x1,y1,'x-','linewidth',1,'markersize',4,'DisplayName','N=10');
hold on;
semilogx(x2,y2,'o-','linewidth',1,'markersize',4,'DisplayName','N=13');
hold on;
semilogx(x3,y3,'o-','linewidth',1,'markersize',4,'DisplayName','N=16');
hold on;
semilogx(x4,y4,'o-','linewidth',1,'markersize',4,'DisplayName','N=20');
hold off;
xlabel ('cyclelength')
ylabel 'number of cycles'
legend('show')