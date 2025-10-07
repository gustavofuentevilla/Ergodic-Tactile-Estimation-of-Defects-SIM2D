close all
clear
clc

%% Create ROS 2 node
node = ros2node("/matlab_pose_publisher", 0);

% Create publisher for PoseStamped messages
pub = ros2publisher(node,...
                    "/target_frame",...
                    "geometry_msgs/PoseStamped");

% Create message structure
msg = ros2message(pub);

%% Trajectory parameters
x_c = (0.15 + 0.44)/2;
y_c = (0.23/2);
radius = 0.1;             % radius of the circle
angular_speed = 0.65;      % radians per second
duration = 10;            % seconds
rate = 100;                % Hz
loop_rate = ros2rate(node, rate);

%% Loop time setup
t = 0;
reset(loop_rate);
while t < duration
    % Update timestamp
    now = ros2time(node, "now");
    msg.header.stamp.sec = int32(now.sec);
    msg.header.stamp.nanosec = uint32(now.nanosec);
    msg.header.frame_id = 'world';  % reference frame

    % Circular trajectory (x = r*cos(θ), y = r*sin(θ))
    theta = angular_speed * t;
    msg.pose.position.x = x_c + radius * cos(theta);
    msg.pose.position.y = y_c - radius * sin(theta);
    msg.pose.position.z = 0.2;

    % % Orientation as yaw (2D), converted to quaternion
    % yaw = theta + pi/2;  % Facing tangent to the circle
    % q = eul2quat([yaw 0 0]);  % [yaw pitch roll] → [qx qy qz qw]
    % 
    % msg.pose.orientation.x = q(2);
    % msg.pose.orientation.y = q(3);
    % msg.pose.orientation.z = q(4);
    % msg.pose.orientation.w = q(1);
    msg.pose.orientation.x = sqrt(2)/2;
    msg.pose.orientation.y = sqrt(2)/2;
    msg.pose.orientation.z = 0.0;
    msg.pose.orientation.w = 0.0;

    % Publish the message
    send(pub, msg);
    disp("published" + t)

    % Wait for the next iteration
    waitfor(loop_rate);
    t = t + 1/rate;
end

%% Prueba trayectoria

tt = (0:0.01:10)';
theta_tmp = 0.65*tt;
x = x_c + radius*cos(theta_tmp);
y = y_c + radius*sin(theta_tmp);

plot(x,y,"--")
axis equal
hold on
comet(x,y)
hold off