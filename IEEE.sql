#create database IEEE14;
use IEEE14;
create table if not exists bus
(
bus int unsigned,
name varchar(20),
loadFlowAreaNo int unsigned,
lossZoneNo int unsigned,
type tinyint unsigned,
finalVolt float(15, 10),
finalAngle float(15, 10),
loadMW float(15, 10),
loadMVAR float(15, 10),
generationMW  float(15, 10),
generationMVAR float(15, 10),
baseKV float(15, 10),
desiredVolt float(15, 10),
maxMVAR float(15, 10),
minMVAR float(15, 10),
G float(15, 10),
B float(15, 10),
remoteCtrlBusNo int unsigned
);

create table if not exists branch
(
tapBusNo int unsigned,
ZbusNo int unsigned,
loadFlowAreaNo int unsigned,
lossZoneNo int unsigned,
circuit int unsigned,
type tinyint unsigned,
R float(15, 10),
X float(15, 10),
B float(15, 10),
lineMVAratingNo1 float(15, 10),
lineMVAratingNo2 float(15, 10),
lineMVAratingNo3 float(15, 10),
ctrlBusNo int unsigned,

side tinyint unsigned,
transformerFinalTurnsRatio float(15, 10),
transformerPhaseShifterFinalAngle float(15, 10),
minTapOrPhaseShift float(15, 10),
maxTapOrPhaseShift float(15, 10),
stepSize float(15, 10),
minVolt float(15, 10),
maxVolt float(15, 10)
);