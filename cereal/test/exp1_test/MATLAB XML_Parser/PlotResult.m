clc
close all
clear all

addpath ../

%% system
stateTrj = parseXML('exp1State.xml',[0 0; 1 0]);
timeTrj  = parseXML('exp1Time.xml',[0 0]);
inputTrj  = parseXML('exp1Input.xml',[0 0]);

figure(1)
plot(stateTrj(1,:), stateTrj(2,:))

figure(2)
subplot(2,1,1)
plot(timeTrj, stateTrj(1,:))
subplot(2,1,2)
plot(timeTrj, stateTrj(2,:))

figure(3)
plot(timeTrj, inputTrj)

%% Sensitivity
stateSensitivityTrj = parseXML('exp1StateSensitivity.xml',[0 0; 1 0]);
timeSensitivityTrj  = parseXML('exp1TimeSensitivity.xml',[0 0]);

figure(4)
subplot(2,1,1)
plot(timeSensitivityTrj, stateSensitivityTrj(1,:))
subplot(2,1,2)
plot(timeSensitivityTrj, stateSensitivityTrj(2,:))
