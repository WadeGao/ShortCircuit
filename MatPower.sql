create database DataSet250k;
use DataSet250k;

create table if not exists branch(
fbus int unsigned default 0 not null,
tbus int unsigned default 0 not null,

r float(15,10),
x float(15,10),
b float(15,10),

rateA float(15,10) default 0,
rateB float(15,10),
rateC float(15,10),

ratio float(15,10),
angle float(15,10) default 0.00,
`status` BIT(1),
angmin float(15,10),
angmax float(15,10),
Pf float(15,10),
Qf float(15,10),
Pt float(15,10),
Qt float(15,10),
mu_Sf float(15,10),
mu_St float(15,10),
mu_angmin float(15,10),
mu_angmax float(15,10),

# primary key (fbus, tbus),
# unique index socket(fbus, tbus),
index socket(fbus, tbus),
index transformer(ratio)
);

create table if not exists generator(
bus int unsigned default 0 not null,
Pg float(15,10),
Qg float(15,10),
Qmax float(15,10),
Qmin float(15,10),
Vg float(15,10),
mBase float(15,10),
`status` BIT(1),
Pmax float(15,10),
Pmin float(15,10),

Pc1 float(15,10) default 0.00,
Pc2 float(15,10) default 0.00,
Qc1min float(15,10) default 0.00,
Qc1max float(15,10) default 0.00,
Qc2min float(15,10) default 0.00,
Qc2max float(15,10) default 0.00,
ramp_agc float(15,10) default 0.00,
ramp_10 float(15,10) default 0.00,
ramp_30 float(15,10) default 0.00,
ramp_q float(15,10) default 0.00,
apf float(15,10) default 0.00,
mu_Pmax float(15,10) default 0.00,
mu_Pmin float(15,10) default 0.00,
mu_Qmax float(15,10) default 0.00,
mu_Qmin float(15,10) default 0.00
);

create table if not exists bus(
bus_i int unsigned not null,
`type` TINYINT unsigned,
Pd float(15,10),
Qd float(15,10),
Gs float(15,10),
Bs float(15,10),
`area` int unsigned default 1,
Vm float(15,10),
Va float(15,10),
baseKV float(15,10),
zone int unsigned default 1,
Vmax float(15,10),
Vmin float(15,10),
lam_P float(15,10),
lam_Q float(15,10),
mu_Vmax float(15,10),
mu_Vmin float(15,10),

primary key (bus_i),
unique index nodeNum(bus_i)
);