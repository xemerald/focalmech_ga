
# This is focalmech_ga's parameter file

# Basic Earthworm setup:
#
MyModuleId          MOD_FMECH_GA    # module id for this instance of shakemap
RingName            HYPO_RING       # shared memory ring for input/output
LogFile             1               # 0 to turn off disk log file; 1 to turn it on
                                    # to log to module log but not stderr/stdout
HeartBeatInterval   15              # seconds between heartbeats

ThreadsNum          4               # number of parallel processing threads, max is 128

RemoveFocalPlot     0               # 0 to keep those plot files; 1 to remove those
                                    # files after posted. Defaults to 0 if this is not setted
MinPickPolarity     3               # min number of pick polarity

# Parameters for main Genetic Algorithm:
#
IterationNum        20              # maximum iteration times, defaults to 20
Population          800             # population number for one generation, defaults to 800
MutateBits          3               # number of bit for mutation once, defaults to 3
ReprodutionRate     0.036           # ratio of population for reproduction, defaults to 0.036 (3.6%)
MutateRate          0.72            # ratio of population for mutation, defaults to 0.72 (72%)

# Directory to create the report files:
#
ReportPath           /home/.../ew/run/focalmechs

# File define the 3D P-wave & S-wave velocity model
#
VelocityModelFile     /home/.../ew/run/params/VPVSMOD.txt

# Post to the other place:
# This function is designed especially for executing external script.
# And it will be called like this: "script_name report_time magnitude trigstations
# result_filename_1"
# If you don't want to use it, please comment it out!
#
# PostScript            /home/.../ew/run/params/post_facebook.py
# PostScriptWithMinMag  /home/.../ew/run/params/post_facebook.py    5.0

# List the message logos to grab from transport ring:
#
#              Installation       Module          Message Types
GetEventsFrom  INST_WILDCARD    MOD_WILDCARD     TYPE_FULL_EVENT
