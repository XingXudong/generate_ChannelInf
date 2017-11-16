function simConsts = simConfig()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 仿真系统配置函数    
% 功能：实现仿真系统参数配置 (上行)   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% 输入：无                               
% 输出：配置参数结构体                    
% 使用示例：simConsts = simConfig()      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%可手动调整部分
simConsts = struct('userPerCell', 10, ...           %每小区用户数              
                   'bandWidth', 10e6, ...           %带宽     Hz
                   'interSiteDist', 500, ...    	%基站间距ISD   m
                   'BS_antenna', 4, ...             %基站发射天线个数 
                   'MS_antenna', 1, ...             %UE发射天线个数                
                   'fc', 2e9, ...                   %载波频率   Hz
                   'Scenario', 'urban_macro', ...   %仿真环境
                   'BsElementPosition', 5, ...      %基站天线间隔，以波长为单位， 一般取 4 - 10
                   'MsElementPosition', 0.5, ...    %移动台天线间隔，以波长维单位
                   'HARQ_num', 0, ...               %重传次数
                   'allocateRbNum', 5, ...          %每次调度分配的RB数
                   'MS_TX_MAX_Power', 23, ...       %移动台最大发射功率 dBm
                   'MsVelocity', 30, ...             %移动台移速  km/h 300/250/200/150/100/50/0     
                   'thermalNoise', -174, ...        %热噪声， dBm/Hz
                   'BsNoiseFigure', 5, ....         %基站接收端噪声系数， dB
                    ...
                    ... %以下是调度的参数
                    'utility', 'PF', ...            %效用函数，maxCI, PF, RR, deltaR
                    'fairFactor', 1, ...            %PF中的公平系数
                    'maxRbGroup', 5, ...            %每个用户能被分配的最大RB Group数
                    ...
                    ... %以下是上行功率控制的参数 
                    'P0', -60, ....                 %参考功率， dBm  
                    'alpha', 0.6, ...               %补偿比
                    ...
                    ... %coordinate set 
                    'setScheme', '1', ...
                    ...
                    ... %以下是信道模型
                    'channelModel', 'ITU'...        %信道模型，ITU模型或3gpp case 1
                 );
 
             
             
simConsts.p=15;                %跟时延大小一样            
             
             
             
%以下部分无须手动设置
simConsts.siteNum = 19;                             %基站个数
simConsts.sectorPerSite = 3;                        %每个基站的扇区（小区）个数

%输入基站数合法性校验
supportSiteNum = [1, 7,19,37,61,91,127,169,217,271,331]; 
if ~find(supportSiteNum == simConsts.siteNum)
    error('wrong site number!');
end
                           
%带宽、码率、FFT点数、RB数的对应关系
table = [                                               ...
           [1.4,  3,    5,   10,   15,   20]*1e6;       ...     %带宽
           [0.5, 1,    2,   4,    6,    8]*3.84e6;      ...     %码率
           128,  256,  512, 1024, 1536, 2048;           ...     %FFT点
           6,    15,   25,  50,   75,   100;            ...     %RB数
           ];  
       
index = table(1,:) == simConsts.bandWidth;

%带宽合法性校验
if(sum(index) == 0 )
    error('wrong band width');
end

chipRate = table(2, index);                        
fftSize = table(3, index);
rbNum = table(4, index);
simConsts.fftSize = fftSize;                        %FFT点数 
simConsts.chipRate = chipRate;                      %码率
simConsts.rbNum  = rbNum;                           %RB(resouce block)数


simConsts.siteWrapNum = simConsts.siteNum * 7;      %环绕后的基站数
simConsts.cellNum = simConsts.siteNum * 3;          %小区数
simConsts.cellWrapNum = simConsts.cellNum * 7;      %环绕后的小区数
simConsts.userNum = simConsts.userPerCell * simConsts.cellNum;  %用户数
simConsts.userWrapNum = simConsts.userNum * 7;                  %环绕后的用户数

simConsts.rbGroupNum = simConsts.rbNum / simConsts.allocateRbNum;