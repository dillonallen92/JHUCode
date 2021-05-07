
function angV=getAngStation(GS,posFixed)
    vec0=[cos(GS.long)*cos(GS.lat);sin(GS.long)*cos(GS.lat);sin(GS.lat)];
    angV=zeros(size(posFixed,1),1);
    
    for i=1:size(posFixed,1)
        angV(i,1)=acos(dot(vec0,posFixed(i,:)')/norm(vec0)/norm(posFixed(i,:)'));
    end

end