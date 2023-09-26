
# This is focalmech_ga's parameter file

# Basic Earthworm setup:
#
MyModuleId          MOD_FMECH_GA    # module id for this instance of shakemap
RingName            HYPO_RING       # shared memory ring for input/output
LogFile             1               # 0 to turn off disk log file; 1 to turn it on
                                    # to log to module log but not stderr/stdout
HeartBeatInterval   15              # seconds between heartbeats

QueueSize           50              # max messages in internal circular msg buffer

RemoveShakeMap      0               # 0 to keep those shakemap files; 1 to remove those
                                    # files after posted. Defaults to 0 if this is not setted

# Directory to create the report files:
#
ReportPath           /home/.../ew/run/focalmechs

# File define the target zone & city boundary in latitude & longtitude
#
3DVelocityModelFile     /home/.../ew/run/params/VPVSMOD.txt

# Post to the other place:
# This function is designed especially for executing external script.
# And it will be called like this: "script_name start_time end_time report_time max_magnitude trigstations
# result_filename_1 [result_filename_2]..."
# If you don't want to use it, please comment it out!
#
# PostScript            /home/.../ew/run/params/post_facebook.py
# PostScriptWithMinMag  /home/.../ew/run/params/post_facebook.py    5.0

# List the message logos to grab from transport ring:
#
#              Installation       Module          Message Types
GetEventsFrom  INST_WILDCARD    MOD_SHAKEMAP     TYPE_FULL_EVENT
