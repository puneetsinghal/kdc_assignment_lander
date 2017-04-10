% criterion to find rotation between frames.
function score=criterion_inertia(p)

global ang_acceleration
global ang_velocity
global NF

score = 0;
inertia_matrix = [p(1) p(4) p(6);p(4) p(2) p(5); p(6) p(5) p(3)];

for i = 1:NF
    torque = inertia_matrix*ang_acceleration(i,:)' + hat(ang_velocity(i,:))*inertia_matrix*ang_velocity(i,:)';
    score = score + torque'*torque;
end

if(sum(p < zeros(1,6)) >= 1)
    score = NaN;
end

end
