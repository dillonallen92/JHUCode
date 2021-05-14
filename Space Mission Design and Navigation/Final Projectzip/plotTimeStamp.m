%plot the state at time
function plotTimeStamp(x,t,index,rE,nB,GS,figNum)

figure(figNum);
hold on;
ax=gca;

%eliminate any old surface plots
indElim=0;
for i=1:length(ax.Children)
    if(strcmp(class(ax.Children(i)),'matlab.graphics.primitive.Surface'))
        indElim=i;
    end
end
if(indElim)
    ax.Children(indElim).delete;
end

img=imread('Blue_Marble_2002.png');
[xu,yu,zu] = sphere(100);
hsurf=surface(xu*rE,yu*rE,zu*rE,'FaceColor','texturemap','EdgeColor','none','Cdata',flipud(img));
direction = [0 0 1];
rotate(hsurf,direction,nB*t(index)*180/pi,[0,0,0]);        %rotate Earth obj through appropriate angle
axis equal;
hs=scatter3(x(index,1),x(index,2),x(index,3),'filled');

%plot the ground station positions
CfixedToIn=fromInToFixed(nB,t(index))';
for i=1:length(GS)
    vec0=rE*[cos(GS(i).long)*cos(GS(i).lat);sin(GS(i).long)*cos(GS(i).lat);sin(GS(i).lat)];
    vecIn=CfixedToIn*vec0;
    scatter3(vecIn(1,1),vecIn(2,1),vecIn(3,1),'filled');
end

end