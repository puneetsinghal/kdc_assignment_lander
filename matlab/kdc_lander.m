clear;
clc;
%% Part a Quaternion calculations using rotation matrix
% rotation matrix is defined by edges of the artifact as coordinate axes
% and edge vectors are calulated using marker locations

data = mrdplot_convert('d00119');

marker_coordinates = data(:,81:104);

artifact_centroid = reshape(mean(reshape(marker_coordinates,[],8),2),10000,[]);

artifact_frames_x = normr(marker_coordinates(:,22:24)-marker_coordinates(:,10:12));
artifact_frames_y = normr(marker_coordinates(:,22:24)-marker_coordinates(:,16:18));
artifact_frames_z = normr(marker_coordinates(:,16:18)-marker_coordinates(:,13:15));

rotation_artifact = [artifact_frames_x, artifact_frames_y, artifact_frames_z];
quaternion_artifact = zeros(size(rotation_artifact,1),4);
% quaternion_artifact(1,:) = [1, 10^-7, 10^-7, 10^-7];
for i = 1: size(rotation_artifact,1)
    R = [artifact_frames_x(i,:)', artifact_frames_y(i,:)', artifact_frames_z(i,:)'];
    quaternion_artifact(i,:) = ...
        qGetQ([0 -1 0; 1 0 0; 0 0 1]*R);
    if (i>1)
        if (any(abs(quaternion_artifact(i,:)-quaternion_artifact(i-1,:)) > 0.2))
            quaternion_artifact(i,:) = -quaternion_artifact(i,:);
        end
    end
end
% figure; plot(quaternion_artifact(1:4000,:))
%% 
% Part a Center of Mass position and velocity
% Used optimization to find center of mass and velocity

% globals
global markers % marker data array NFx(NM*3)
global NF % number of frames

markers = [linspace(0,1000,11)', marker_coordinates(linspace(1,1001,11),:)];

NF = 10; % number of frames (including initial frame)

% set options for fminunc()
% options = optimset();
options = optimset('MaxFunEvals',1000000);

% p0 is the intitial parameter vector
n_p = 6; % com 3, v 3
p0 = zeros(1,n_p);

% do optimization
[com, ~, ~] = fminunc(@criterion,p0,options);
marker1_offset = marker_coordinates(1,1:3)'-com(1:3)';
com_pos = com(1:3);
com_vel = com(4:end);
com = zeros(size(marker_coordinates,1),3);
com(1,:) = com_pos;
for i = 2:size(marker_coordinates,1)
    %     com(i,:) = marker_coordinates(i,1:3) - ([0 -1 0; 1 0 0; 0 0 1]*reshape(rotation_artifact(i,:),3,3)*marker1_offset)';
    com(i,:) = com_pos + i*1/100*(com_vel);
end

% Print data in Problem_2_0.dat file
fid = fopen('../problem_2_0.dat','wt');
fprintf(fid, 'cx\tcy\tcz\tq_scalar\tq_x\tq_y\tq_z\n');
for i = 1:1001
    fprintf(fid,'%f\t%f\t%f\t%f\t%f\t%f\t%f\n',com(i,1),com(i,2),com(i,3),...
        quaternion_artifact(i,1),quaternion_artifact(i,2),...
        quaternion_artifact(i,3),quaternion_artifact(i,4));
end

%% Part b
% Angular velocity at all times. Calculated using the equation: R_dot =
% R*(w_matrix)

angular_velocity = zeros(size(marker_coordinates,1)-1,3);
dt = 1/100;
for i = 1 : size(marker_coordinates,1)-1
    rot_dot = (rotation_artifact(i+1,:) - rotation_artifact(i,:))/dt;
    ang_matrix = (reshape(rotation_artifact(i,:),3,3))\reshape(rot_dot,3,3);
    angular_velocity(i,:)= unhat(ang_matrix);
end

% Print data in Problem_2_1.dat file
fid = fopen('../problem_2_1.dat','wt');
fprintf(fid, 'wx\twy\twz\n');
for i = 1:1001
    fprintf(fid,'%f\t%f\t%f\n',angular_velocity(i,1), angular_velocity(i,2), angular_velocity(i,3));
end

%% Part c
% Angular acceleration calculations using finite differences

angular_acceleration = zeros(size(marker_coordinates,1)-2,3);
dt = 1/100;
for i = 1:size(marker_coordinates,1)-2
    angular_acceleration(i,:) = (angular_velocity(i+1,:)-angular_velocity(i,:))/dt;
end

% Print data in Problem_2_2.dat file
fid = fopen('../problem_2_2.dat','wt');
fprintf(fid, 'ax\tay\taz\n');
for i = 1:1001
    fprintf(fid,'%f\t%f\t%f\n',angular_acceleration(i,1), angular_acceleration(i,2), angular_acceleration(i,3));
end

%% Part d
% Inertia calculations. Finding Null matrix using SVD decompositions
NF = 1000;
ang_acceleration = [smooth(angular_acceleration(:,1)),...
    smooth(angular_acceleration(:,2)), smooth(angular_acceleration(:,3))];

A = [];
for j = 10:NF
    i = j; %randi(NF);
    temp = [ang_acceleration(i,1), 0, 0, ang_acceleration(i,2), ang_acceleration(i,3), 0;...
        0, ang_acceleration(i,2), 0, angular_acceleration(i,1), 0, angular_acceleration(i,3);...
        0, 0, ang_acceleration(i,3), 0, ang_acceleration(i,1), ang_acceleration(i,2)] +...
        [0, -angular_velocity(i,3) angular_velocity(i,2);...
        angular_velocity(i,3), 0, -angular_velocity(i,1);...
        -angular_velocity(i,2), angular_velocity(i,1), 0]*...
        [angular_velocity(i,1), 0, 0, angular_velocity(i,2), 0, angular_velocity(i,3);...
        0, angular_velocity(i,2), 0, angular_velocity(i,1), 0, angular_velocity(i,3);...
        0, 0, angular_velocity(i,3), 0, angular_velocity(i,1), angular_velocity(i,2)];
    A = [A; temp];
end
[U,S,V] = svd(A);
Inertia = [V(1,end) V(4,end) V(6,end);V(4,end) V(2,end) V(5,end); V(6,end) V(5,end) V(3,end)];
display(Inertia);

%% Part E: Future trajectory

dt = 1/100;
ang_vel_previous = angular_velocity(1001,:)';
com_future = zeros(1000,3);
quat = zeros(1000,4);
R = quat2rotm(quaternion_artifact(1001,:));
fid = fopen('../problem_2_3.dat','wt');
fprintf(fid, 'cx\tcy\tcz\tq_scalar\tq_x\tq_y\tq_z\n');

for i =1:1000
    com_future(i,:) = com_pos+com_vel*(i*dt+10);
    ang_vel_next = -Inertia\(cross(ang_vel_previous, Inertia*ang_vel_previous))*dt + ang_vel_previous;
    dR = dt*R*hat(ang_vel_previous);
    R = R + dR;
    [U,S,V]=svd(R);
    rotation_future{i} = U*V';
    quat(i,:) = rotm2quat(R);
    fprintf(fid,'%f\t%f\t%f\t%f\t%f\t%f\t%f\n',com_future(i,1),com_future(i,2),com_future(i,3),...
        quat(i,1),quat(i,2),...
        quat(i,3),quat(i,4));
    ang_vel_previous = ang_vel_next;
end
fclose(fid);

%% Part 3
% calculating and saving data for landing using original com position and
% quaternion for first 10 seconds and then used the same variables from
% future trajectory (Part 2 E)

NF = 2000;
com_world = zeros(NF, 3);
q = zeros(NF,4);
for i =1:1000
    R = [artifact_frames_x(i,:)', artifact_frames_y(i,:)', artifact_frames_z(i,:)'];
    q(i,:) = ...
        qGetQ([0 -1 0; 1 0 0; 0 0 1]*R*[0 -1 0; 1 0 0; 0 0 1]);
    if (i>1)
        if (any(abs(q(i,:)-q(i-1,:)) > 0.2))
            q(i,:) = -q(i,:);
        end
    end
    com_world_temp = ([0 -1 0 0; 1 0 0 11; 0 0 1 4; 0 0 0 1]*[com(i,:), 1]')';
    %     com_world_temp = ([1 0 0 0; 0 1 0 11; 0 0 1 4; 0 0 0 1]*[com(i,:), 1]')';
    
    com_world(i,:) = com_world_temp(1:3);
    if (i<800)
        offset_assignment = [0 3.5 0]';
    else
        offset_assignment = [0 3.3 -0.5]';
    end
    offset_world = [0 -1 0; 1 0 0; 0 0 1]*reshape(rotation_artifact(i,:),3,3)*offset_assignment;
    com_world(i,:) = com_world(i,:) + offset_world';
end

for i =1001:2000
    
    q(i,:) = ...
        qGetQ(rotation_future{i-1000}*[0 -1 0; 1 0 0; 0 0 1]);
    if (i>1)
        if (any(abs(q(i,:)-q(i-1,:)) > 0.2))
            q(i,:) = -q(i,:);
        end
    end
    com_world_temp = ([0 -1 0 0; 1 0 0 11; 0 0 1 4; 0 0 0 1]*[com_future(i-1000,:), 1]')';
    %     com_world_temp = ([1 0 0 0; 0 1 0 11; 0 0 1 4; 0 0 0 1]*[com(i,:), 1]')';
    
    com_world(i,:) = com_world_temp(1:3);
    if (i<800)
        offset_assignment = [0 3.2 0]';
    else
        offset_assignment = [0.05 3.2 0]';
    end
    offset_world = rotation_future{i-1000}*offset_assignment;
    com_world(i,:) = com_world(i,:) + offset_world';
end


fid = fopen('../lander/q_scalar.txt','wt');
fprintf(fid,'q_scalar\n');
for i =1: NF
    fprintf(fid, '%f\n', q(i,1));
end
fclose(fid);

fid = fopen('../lander/q_x.txt','wt');
fprintf(fid,'q_x\n');
for i =1: NF
    fprintf(fid, '%f\n', q(i,2));
end
fclose(fid);

fid = fopen('../lander/q_y.txt','wt');
fprintf(fid,'q_y\n');
for i =1: NF
    fprintf(fid, '%f\n', q(i,3));
end
fclose(fid);

fid = fopen('../lander/q_z.txt','wt');
fprintf(fid,'q_z\n');
for i =1: NF
    fprintf(fid, '%f\n', q(i,4));
end
fclose(fid);


fid = fopen('../lander/COM_X.txt','wt');
fprintf(fid,'COM_X\n');
for i =1: NF
    fprintf(fid, '%f\n',com_world(i,1));
end
fclose(fid);

fid = fopen('../lander/COM_Y.txt','wt');
fprintf(fid,'COM_Y\n');
for i =1: NF
    fprintf(fid, '%f\n',com_world(i,2));
end
fclose(fid);

fid = fopen('../lander/COM_Z.txt','wt');
fprintf(fid,'COM_Z\n');
for i =1: NF
    fprintf(fid, '%f\n',com_world(i,3));
end
fclose(fid);

fid = fopen('../lander/COM_Lander.dat','wt');
fprintf(fid,'COM\n');
for i =1: NF
    fprintf(fid, '%f\t%f\t%f\n',com(i,1), com(i,2), com(i,3));
end
fclose(fid);


fid = fopen('../lander/ang_X.txt','wt');
fprintf(fid,'ang_X\n');
for i =1: NF
    fprintf(fid, '%f\n',angular_velocity(i,1));
end
fclose(fid);

fid = fopen('../lander/ang_Y.txt','wt');
fprintf(fid,'ang_Y\n');
for i =1: NF
    fprintf(fid, '%f\n',angular_velocity(i,2));
end
fclose(fid);

fid = fopen('../lander/ang_Z.txt','wt');
fprintf(fid,'ang_Z\n');
for i =1: NF
    fprintf(fid, '%f\n',angular_velocity(i,3));
end
fclose(fid);