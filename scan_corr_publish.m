[n, m] = size(out.res);
opt.fontsize = 18;

% prepare a figure
w1 = 0.1;
w2 = 0.57;
s = 0.05;
xpos(1) = 0.13;
xpos(2) = xpos(1)+w1+s;

h1 = 0.1;
h2 = 0.6;
ypos(1) = 0.1;
ypos(2) = xpos(1)+h1+s;

col1 = [     0    0.4470    0.7410];
col2 = [0.8500    0.3250    0.0980];

figpos = [50 50 830 800];
figure('Position',figpos);
hdl_fig = gcf;

% set the paper position to 1 inch per 100 pixels
set(hdl_fig,'PaperUnits','inches','PaperPosition',figpos.*[0 0 1e-2 1e-2])
set(hdl_fig,'PaperSize',figpos(3:4).*[1e-2 1e-2])

hdl_a(1) = axes('Position',[xpos(2),ypos(2),w2,h2],'FontSize',opt.fontsize);
hdl_i = imagesc(out.res);
set(hdl_a(1),'Box','On','DataAspectRatio',[1 1 1]);
set(hdl_a(1),'Xtick',{},'Ytick',{});
colorbar
% colormap(cmap);
set(hdl_a(1),'Position',[xpos(2),ypos(2),w2,h2]);
% plot the ROI
hdl_roi = patch('Vertices',nan(4,2),'Faces',1:4,'FaceColor','none','EdgeColor',col1,'Parent',hdl_a(1));

hdl_a(2) = axes('Position',[xpos(1),ypos(2),w1,h2],'NextPlot','Add','Box','On','FontSize',opt.fontsize,'YDir','reverse');
hdl_p(1) = plot(nan,nan,'-','Color',col1);
hdl_p(2) = plot(nan,nan,'-','Color',col2);
legend(hdl_p(1:2),{'g_{hx}','g_{hy}'},'Orientation','horizontal','Location','NorthOutside','FontSize',opt.fontsize);
set(hdl_a(2),'Position',[xpos(1),ypos(2),w1,h2]);
xlabel('$u_h$ [px]','Interpreter','Latex','FontSize',opt.fontsize)
ylabel('$y$ [px]','Interpreter','Latex','FontSize',opt.fontsize)

hdl_a(3) = axes('Position',[xpos(2),ypos(1),w2,h1],'NextPlot','Add','Box','On','FontSize',opt.fontsize);
hdl_p(3) = plot(nan,nan,'-','Color',col2);
hdl_p(4) = plot(nan,nan,'-','Color',col1);
legend(hdl_p(3:4),{'g_{vx}','g_{vy}'},'Orientation','vertical','Location','EastOutside','FontSize',opt.fontsize);
ylabel('$u_v$ [px]','Interpreter','Latex','FontSize',opt.fontsize)
xlabel('$x$ [px]','Interpreter','Latex','FontSize',opt.fontsize)
set(hdl_a(3),'Position',[xpos(2),ypos(1),w2,h1]);


set(hdl_roi,'Vertices',nan(4,2));
set(hdl_a(2),'YLim',[1 n]);
set(hdl_a(3),'XLim',[1 m]);
set(hdl_p(1),'YData',1:n,'XData',out.uh(:,1));
set(hdl_p(2),'YData',1:n,'XData',out.uh(:,2));
set(hdl_p(3),'XData',1:m,'YData',out.uv(1,:));
set(hdl_p(4),'XData',1:m,'YData',out.uv(2,:));


