function plotOrb(x,figNum)

figure(figNum);
hold on;
hp=plot3(x(:,1),x(:,2),x(:,3));
axis equal;
hp.LineWidth=3;
ax=gca;
ax.FontSize=20;
ax.FontName='Times New Roman';

ax.XLabel.String='x (km)';
ax.YLabel.String='y (km)';
ax.ZLabel.String='z (km)';

ax.XLabel.Interpreter='latex';
ax.YLabel.Interpreter='latex';
ax.ZLabel.Interpreter='latex';
box on;

end