function [altV, latV, longV] = posFixedToAltLatLong(posFixed, rE)
    % Implementation of posFixedToAltLatLong based on page 39
    % of TSB book.
    for i = 1 : length(posFixed)
       x = posFixed(i,1);
       y = posFixed(i,2);
       z = posFixed(i,3);
       r = sqrt(x^2 + y^2 + z^2);
       longV(i) = atan2(y,x);
       latV(i) = atan2(z,sqrt(x^2+y^2));
       altV = r - rE;
    end
    latV = latV';
    longV = longV';
    altV = altV';
end