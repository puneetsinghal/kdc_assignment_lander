function Skew_Mat = hat(p)


p1 = p(1);
p2 = p(2);
p3 = p(3);

Skew_Mat = [0 -p3 p2 ; p3 0 -p1 ; -p2 p1 0];


end

