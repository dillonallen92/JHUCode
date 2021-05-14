
function angV=getAngStation(GS,posFixed)
    vec0=[cos(GS.long)*cos(GS.lat);sin(GS.long)*cos(GS.lat);sin(GS.lat)];
    angV=zeros(size(posFixed,1),1);
    
    
    
    
    for i=1:size(posFixed,1)
        rsat_station=posFixed(i,:)'-vec0;
       % angV(i,1)=acos(dot(vec0,posFixed(i,:)')/norm(vec0)/norm(posFixed(i,:)'));
        
         angV(i,1)=acos(dot(vec0,rsat_station)/norm(vec0)/norm(rsat_station));
    end

end