% Use optimization to find rotation between frames.

% globals
global markers % marker data array NFx(NM*3)
global NF % number of frames

load markers

NF = 10; % number of frames (including initial frame)

% set options for fminunc()
% options = optimset();
options = optimset('MaxFunEvals',1000000);

% p0 is the intitial parameter vector
n_p = 4*(NF-1); % (NF-1) * (axis 3 + angle 1)
p0(n_p) = 0;
for i = 1:n_p
 p0(i) = 0;
end
% set up initial axis
index = 1;
for i = 1:NF-1;
 p0(index) = 1;
 index = index + 4;
end

% do optimization
[answer,fval,exitflag]=fminunc(@criterion_omega,p0,options)

