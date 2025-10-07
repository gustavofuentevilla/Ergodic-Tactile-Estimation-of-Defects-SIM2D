close all
clear
clc

%% Define node

nodeMATLAB = ros2node("nodeMATLAB", 0);

%% Check nodes

ros2 node list

%% Check topics and message types

ros2 topic list -t

%% Explore message data structure

ros2 msg show geometry_msgs/Twist
ros2 msg show geometry_msgs/Vector3
ros2 msg show turtlesim/Pose %% Error because the msg type is not recognized by matlab, for that we need custom msgs

%% List all available message types 

ros2 msg list

%% Define Publisher

publisher = ros2publisher(nodeMATLAB, ...
                          "/turtle1/cmd_vel");   %topic

%% Create ROS2 Message and assign values

msg = ros2message(publisher);
msg.linear.x = 2; 
msg.angular.z = 1;

%% Specify execution rate for the node

ratePublisher = ros2rate(nodeMATLAB, 1); %1 Hz

%% Loop

reset(ratePublisher)

for i = 1:500
    send(publisher, msg); % Publish the message
    disp(i + "Published: " + msg.data)
    waitfor(ratePublisher); % Wait for the next iteration
end