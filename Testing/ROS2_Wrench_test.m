close all
clear
clc

%% Create ROS 2 node
nodeWrench = ros2node("/matlab_wrench_publisher", 0);

%% Create publisher for PoseStamped messages
pub = ros2publisher(nodeWrench,...
                    "/cartesian_compliance_controller/target_wrench",...
                    "geometry_msgs/WrenchStamped");

% Create message structure
msg = ros2message(pub);

%% parameters
duration = 60;            % seconds
rate = 100;                % Hz
loop_rate = ros2rate(nodeWrench, rate);

%% Loop time setup
t = 0;
reset(loop_rate);
while t < duration
    % Update timestamp
    now = ros2time(nodeWrench, "now");
    msg.header.stamp.sec = int32(now.sec);
    msg.header.stamp.nanosec = uint32(now.nanosec);
    msg.header.frame_id = 'ur5e_tool0';  % reference frame

    % Wrench
    msg.wrench.force.x = 0.0;
    msg.wrench.force.y = 0.0;
    msg.wrench.force.z = 5.0;
    msg.wrench.torque.x = 0.0;
    msg.wrench.torque.y = 0.0;
    msg.wrench.torque.z = 0.0;

    % Publish the message
    send(pub, msg);
    disp("published " + t)

    % Wait for the next iteration
    waitfor(loop_rate);
    t = t + 1/rate;
end