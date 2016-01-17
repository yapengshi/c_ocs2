clc
close all
clear all

addpath ../

%%
h = tf(1,[1 2 1]);
step(h,10)

%% Model Subsystem
stateTrj = parseXML('secondOrderState.xml',[0 0; 1 0]);
timeTrj  = parseXML('secondOrderTime.xml',[0 0]);

hold on
plot(timeTrj, stateTrj(2,:), '*r')
hold off
