function q_out = rot2quat_without_gimbal(r, q_prev)
 
% q_alien[4]={q_prev(1),q_prev(2),q_prev(3),q_prev(4)};
 scalar1 = 0.5*sqrt(r(1,1)+r(2, 2)+r(3, 3)+1);
 scalar2 = -scalar1;
 eX1 = real(0.5 * sign(r(3, 2)-r(2, 3))*sqrt(r(1,1)-r(2, 2)-r(3, 3)+1));
 eX2 = -eX1;
 eY1 = real(0.5*sign(r(1, 3)-r(3, 1))*sqrt(r(2, 2)-r(3, 3)-r(1,1)+1));
 eY2 = -eY1;
 eZ1 = real(0.5*sign(r(2, 1)-r(1,2))*sqrt(r(3, 3)-r(1,1)-r(2, 2)+1));
 eZ2 = -eZ1;

 q_out1 = [scalar1, eX1, eY1, eZ1];
 q_out2 = [scalar2, eX2, eY2, eZ2];


if (sign(q_prev(1)) == sign(q_out1(1)))
    q_out(1)=q_out1(1);
else
    q_out(1)=q_out2(1);
end
if (sign(q_prev(2)) == sign(q_out1(2)))
    q_out(2) = q_out1(2);
else
    q_out(2)=q_out2(2);
end
if (sign(q_prev(3)) == sign(q_out1(3)))
    q_out(3) = q_out1(3);
else
    q_out(3) = q_out2(3);
end
if (sign(q_prev(4)) == sign(q_out1(4)))
    q_out(4)=q_out1(4);
else
    q_out(4)=q_out2(4);
end
end