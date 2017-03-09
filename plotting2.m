
h1 = scatter(x1,y1,'s','displayName','k=0, \alpha_{1}=0.202')
hold on
h2 = scatter(x2,y2,'d','displayName','k=0.5 \alpha_{1}=0.149')
hold on
a = 0.2024; b =  -0.102; x = 1:16;
h3 = plot(x, a*x+b);
hold on
a = 0.1422; b =   -0.133; x = 1:16;
h4 = plot(x, a*x+b);
hold off
xlabel 'N'
ylabel 'log<N_{1}>'
h_legend=legend([h1 h2])
set(h_legend,'FontSize',14);