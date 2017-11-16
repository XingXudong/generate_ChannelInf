function [H, delayChip] = generateFastFadingChannel() 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 生成快衰信道函数    
% 功能：使用3GPP TR25.996记录的SCM信道模型和SCM代码,生成快衰信道函数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 输入:  simConsts, 仿真配置结构体
% 输出：H delayChip
%      1. H, 快衰信道，[cellNum, Nr, Nt, 6, userNum]维
%      2. delayChip, 各径延迟[cellNum, userNum, 6]维，以码片为单位
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

simConsts = simConfig();

%添加当前路径下的所有子目录
path = 'C:\Users\WTI\Desktop\信道产生及预测模块';
addpath(genpath(path));

%参数初始化
NumTimeSamples = 350;                                         %每条链路采样点数
MsVelocity = simConsts.MsVelocity / 3.6;                    %移动台速率, m/s
chipRate = simConsts.chipRate;                              %码率/采样率
userNum = simConsts.userNum;                                %用户数
cellNum = simConsts.cellNum;                                %小区数
NumMsElements = simConsts.MS_antenna;                       %移动台天线数
NumBsElements = simConsts.BS_antenna;                       %基站天线数
bandWidth = simConsts.bandWidth;

%参数计算      
DelaySamplingInterval =   1 / (chipRate * 16);        %时延采样间隔/分辨率

%SCM参数设置
scmpar = scmparset;                                         %基本设置，默认不包括路径损耗和阴影衰落
scmpar.SampleDensity = 9;                                   %采样密度，单位为sample/half wavelength
% scmpar.UniformTimeSampling = 'yes';
scmpar.NumTimeSamples = NumTimeSamples;                     %采样点数
scmpar.DelaySamplingInterval = DelaySamplingInterval;       %采样间隔
scmpar.NumBsElements = NumBsElements;                       %基站天线数量
scmpar.NumMsElements = NumMsElements;                       %移动台天线数量
scmpar.CenterFrequency = simConsts.fc;                      %载波频率
scmpar.Scenario = simConsts.Scenario;                       %场景
scmpar.BsUrbanMacroAS = 'fifteen';                          %BS 平均角度扩展 
%链路设置
linkpar = linkparset(userNum);                              %基本设置
linkpar.MsVelocity = MsVelocity*ones(1, userNum);           %移动台速率, m/s
%天线参数设置
antpar = antparset;                                         %天线基本设置
antpar.BsElementPosition = simConsts.BsElementPosition;     %基站天线间隔
antpar.MsElementPosition = simConsts.MsElementPosition;     %移动台天线间隔

%生成信道信息
H = zeros(NumTimeSamples, NumMsElements, NumBsElements, 6, userNum);
delayChip = zeros(NumTimeSamples, userNum, 6);
[H1, delayChip1,  ~] =  scm(scmpar, linkpar, antpar);   %生成该小区到所有用户的SCM信道
for i = 1:NumTimeSamples    %NumTimeSamples/10为子帧数
    %     H1 = squeeze(mean(H1, 4));                              %所有采样点取平均
    delayChip1 = round(delayChip1 * chipRate);              %快衰落信道中的多径延迟，并转化为码片数
    H(i, :, :, :, :) = H1(:,:,:,i,:);
    delayChip(i, :, :) = delayChip1;
end
 
%将信道信息临时存储一下
saveDir = 'G:\maozezhong\H.mat';
save(saveDir,'H');






    
    
 
