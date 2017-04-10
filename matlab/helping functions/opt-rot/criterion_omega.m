% criterion to find rotation between frames.
function score=criterion_omega(p)

global markers
global NF

NM = 8; % number of markers
dt = 0.01;

% initialize arrays
axes(NF,3) = 0;
sa(NF) = 0;
ca(NF) = 0;

% pull out parameters
com = [ -11 0 0 ];
v = [ 0.02 -0.05 0.01 ];

j = 1;
for i = 2:NF
 axes(i,1) = p(j);
 axes(i,2) = p(j+1);
 axes(i,3) = p(j+2);
 sa(i) = sin( p(j+3) );
 ca(i) = cos( p(j+3) );
 j = j + 4;
end

score = 0;

% normalize axis
for i = 2:NF
 axis =  axes(i,:);
 norm = axis * axis';
% check for zero w, avoid divide by zero
 if norm < 1e-7
  norm = 1e-7;
 end
 axis = axis/sqrt(norm);
% add a pinch of penalty on w to normalize it to length 1.
 score = score + 0.001*(norm - 1)*(norm - 1);
 axes(i,1) = axis(1);
 axes(i,2) = axis(2);
 axes(i,3) = axis(3);
end

% Use Rodrigues' rotation formula
% http://en.wikipedia.org/wiki/Axis-angle_representation
% v*ca + axisxv*sa + axis*(axis*v')(1-ca)

for i = 2:NF
 for j = 1:NM
  % last +1 to skip initial count variable
  m1 = markers(1,(3*(j-1)+1+1):(3*(j-1)+3+1)) - com;
  mi = markers(i,(3*(j-1)+1+1):(3*(j-1)+3+1));
  mi_est = m1*ca(i) + cross(axes(i,1:3),m1)*sa(i) + axes(i,1:3)*(axes(i,1:3) * m1')*(1-ca(i)) + com + v*dt*markers(i,1);
%  (mi - mi_est)
  score = score + (mi - mi_est)*(mi - mi_est)';
 end
end

end
