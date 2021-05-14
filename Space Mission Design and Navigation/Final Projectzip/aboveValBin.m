
function binV=aboveValBin(paramV,val0)
    
    binV=zeros(length(paramV),1);
    for i=1:length(paramV)
       if(paramV(i,1)>val0)
            binV(i,1)=1;
       end
    end

end