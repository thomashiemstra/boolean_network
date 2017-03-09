sz = 100;
scatter(x1,y1,sz,'x','displayName','\eta=0.4')
hold on
scatter(x2,y2,sz,'+','displayName','\eta=0.6')
hold off
xlabel 'N'
ylabel '<l^{\prime}>'
h_legend = legend('show')
set(h_legend,'FontSize',14);