.ALIASES
C_C1            C1(1=VSEI 2=VBATT ) CN @BATTERY_2_STATES_CURRENT_INPUT.SCHEMATIC1(sch_1):INS737@ANALOG.C.Normal(chips)
V_V1            V1(+=VOCV -=0 ) CN @BATTERY_2_STATES_CURRENT_INPUT.SCHEMATIC1(sch_1):INS883@SOURCE.VDC.Normal(chips)
L_L1            L1(1=VL0 2=VSEI ) CN @BATTERY_2_STATES_CURRENT_INPUT.SCHEMATIC1(sch_1):INS769@ANALOG.L.Normal(chips)
R_R1            R1(1=VOCV 2=VL0 ) CN @BATTERY_2_STATES_CURRENT_INPUT.SCHEMATIC1(sch_1):INS791@ANALOG.R.Normal(chips)
R_R2            R2(1=VSEI 2=VBATT ) CN @BATTERY_2_STATES_CURRENT_INPUT.SCHEMATIC1(sch_1):INS851@ANALOG.R.Normal(chips)
I_I1            I1(+=VBATT -=GND ) CN @BATTERY_2_STATES_CURRENT_INPUT.SCHEMATIC1(sch_1):INS1066@SOURCE.IPULSE.Normal(chips)
R_R3            R3(1=0 2=GND ) CN @BATTERY_2_STATES_CURRENT_INPUT.SCHEMATIC1(sch_1):INS1205@ANALOG.R.Normal(chips)
_    _(gnd=GND)
_    _(Vbatt=VBATT)
_    _(VL0=VL0)
_    _(Vocv=VOCV)
_    _(Vsei=VSEI)
.ENDALIASES
