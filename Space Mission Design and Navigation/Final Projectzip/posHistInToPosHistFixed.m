function posFixed=posHistInToPosHistFixed(posIn,t,nB)
posFixed=zeros(size(posIn));                        %initialize
for i=1:length(t)
    C=fromInToFixed(nB,t(i));
    posFixedTemp=C*posIn(i,:)';
    posFixed(i,:)=posFixedTemp';
end
end