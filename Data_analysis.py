import csv
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
#import timeit
#import torch

path = "C:/Users/qfa20/OneDrive - University of Cambridge/Papers/Opti space depth thermal conductivity/CalculationsV5/"

def lecture(filename):
    a = []
    data = []
    with open(path + filename,'r') as f:
        data = []
        for x in f:
            data.append(x.strip('\n'))
        f.close()
    data_output = [[] for _ in range(len(data))]
    for i in range(len(data)):
        data[i] = data[i].split(" ")
        for j in range(len(data[i])):
            if data[i][j] != '':
                data_output[i].append(float(data[i][j]))
    return(data_output)

#to launch the function split, call a = lecture()

list_TL_min_summary = []
list_TL_max_summary = []
list_TR_min_summary = []
list_TR_max_summary = []

list_HL_min_summary = []
list_HL_max_summary = []
list_HR_min_summary = []
list_HR_max_summary = []

list_thermalcond_summary = []
list_Rspacing_summary = []
list_Rdepth_depth = []
list_Rpower_summary = []

number_of_files = 336

for i in range(number_of_files):
    a = []
    a = lecture("Tdepth." + str(i) + ".dat")
    list_TL_min = []
    list_TL_max = []
    list_TR_min = []
    list_TR_max = []
    list_HL_min = []
    list_HL_max = []
    list_HR_min = []
    list_HR_max = []
    list_thermalcond = []
    list_Rspacing = []
    list_Rdepth = []
    list_Rpower = []
    for j in range(1,len(a)):
        list_TL_min.append(a[j][1])
        list_TL_max.append(a[j][1])
        list_TR_min.append(a[j][2])
        list_TR_max.append(a[j][2])
        list_HL_min.append(a[j][3])
        list_HL_max.append(a[j][3])
        list_HR_min.append(a[j][4])
        list_HR_max.append(a[j][4])
        list_thermalcond.append(a[j][5])
        list_Rspacing.append(a[j][6])
        list_Rdepth.append(a[j][7])
        list_Rpower.append(a[j][8])
    list_TL_min_summary.append(min(list_TL_min))
    list_TL_max_summary.append(max(list_TL_max))
    list_TR_min_summary.append(min(list_TR_min))
    list_TR_max_summary.append(max(list_TR_max))
    list_HL_min_summary.append(min(list_HL_min))
    list_HL_max_summary.append(max(list_HL_max))
    list_HR_min_summary.append(min(list_HR_min))
    list_HR_max_summary.append(max(list_HR_max))
    list_thermalcond_summary.append(max(list_thermalcond))
    list_Rspacing_summary.append(max(list_Rspacing))
    list_Rdepth_depth.append(max(list_Rdepth))
    list_Rpower_summary.append(max(list_Rpower)/(10.0*3600.0))
    
for i in range(number_of_files):
    list_TL_min_summary[i] = round(list_TL_min_summary[i],2)
    list_TL_max_summary[i] = round(list_TL_max_summary[i],2)
    list_TR_min_summary[i] = round(list_TR_min_summary[i],2)
    list_TR_max_summary[i] = round(list_TR_max_summary[i],2)
    list_HL_min_summary[i] = round(list_HL_min_summary[i],2)
    list_HL_max_summary[i] = round(list_HL_max_summary[i],2)
    list_HR_min_summary[i] = round(list_HR_min_summary[i],2)
    list_HR_max_summary[i] = round(list_HR_max_summary[i],2)
    list_thermalcond_summary[i] = round(list_thermalcond_summary[i],2)
    list_Rspacing_summary[i] = round(list_Rspacing_summary[i],2)
    list_Rdepth_depth[i] = round(list_Rdepth_depth[i],3)
    list_Rpower_summary[i] = round(list_Rpower_summary[i],2)
    
outputfile = open(path + 'Summary.dat', "w")
outputfile.write("Case [--]" + "\t" + "TL min [C]" + "\t" + "TL max [C]" + "\t" + "TR min [C]" + "\t" + "TR max [C]" + "\t" + "HL min [W/m**2]" + "\t" + "HL max [W/m**2]" + "\t" + "HR min [W/m**2]" + "\t" + "HR max [W/m**2]" + "\t" + "Thermal cond [W/m/K]" + "\t" + "Ribbons spacing [m]" + "\t" + "Ribbons depth [m]" + "\t" + "Ribbons activation duration [h]" + "\n")
for i in range(number_of_files):
    outputfile.write(str(i) + "\t" + str(list_TL_min_summary[i]) + "\t" + str(list_TL_max_summary[i]) + "\t" + str(list_TR_min_summary[i]) + "\t" + str(list_TR_max_summary[i]) + "\t" + str(list_HL_min_summary[i]) + "\t" + str(list_HL_max_summary[i]) + "\t" + str(list_HR_min_summary[i]) + "\t" + str(list_HR_max_summary[i]) + "\t" + str(list_thermalcond_summary[i]) + "\t" + str(list_Rspacing_summary[i]) + "\t" + str(list_Rdepth_depth[i]) + "\t" + str(list_Rpower_summary[i]) + "\n")
outputfile.close()

