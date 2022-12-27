import os
import math
from datetime import datetime

###########################################
# User setting
###########################################

InputName="myMVTXhits.txt"

OutputName="Hitdata_0Tesla_eta0.txt"

NumberOfZPosition = 20

NumberOfAngle = 180

TrackNumOffset = 0

###########################################
EventsPerEta = NumberOfZPosition * NumberOfAngle

fresult=open(OutputName,'w')

fs= open(InputName,'r')
linedat=fs.readlines()
fresult.write("eta, phi, event num., track num., layer ID, stave ID, chip ID, row, column, global_x, global_y, global_z, trackType, vertex_x, vertex_y, vertex_z\n")
for headline1 in range(0,len(linedat)):
    if headline1 == 0 :
        continue
    eventnum = int(linedat[headline1].split(sep=',')[0])
    if eventnum < EventsPerEta*1 :
        fresult.write("-0.002, ")
    elif eventnum < EventsPerEta*2 :
        fresult.write("-0.001, ")
    elif eventnum < EventsPerEta*3 :
        fresult.write("0.000, ")
    elif eventnum < EventsPerEta*4 :
        fresult.write("0.001, ")
    else:
        fresult.write("0.002, ")
    tempalign=divmod(divmod(eventnum,NumberOfZPosition)[0],NumberOfAngle)
    fresult.write(str(tempalign[1]*2-178)+", "+str(int(linedat[headline1].split(',',maxsplit=1)[0])+TrackNumOffset).zfill(6)+linedat[headline1].split(',',maxsplit=1)[1])
fs.close()
fresult.close()

