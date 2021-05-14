function binV=aboveBelowBin(paramV,valLb,valUb)
    binV=zeros(length(paramV),1);
    for i=1:length(paramV)
       if(paramV(i,1)>valLb&&paramV(i,1)<valUb)
            binV(i,1)=1;
       end
    end
end