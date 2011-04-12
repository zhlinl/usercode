#!/usr/bin/python

import sys
import os
import os.path
import time

"""
runs the systemtics fits
"""

files = list()

print sys.argv
if len(sys.argv) != 3 : 
    print 'enter two arguments' \
        '\n\tfirst argument 1-3 to select 1s/2s/3s'  \
        '\n\tsecond argument 0-2 to select ana mode: ptdiff/1ybin, ptdiff/2ybin, ydiff'
    sys.exit(10) 

if not os.path.isfile("oniafitter.C") :
    print """fitting program is missing in directory
          ln -s ../oniafitter.C ."""
    sys.exit(10) 


pol  = 1
peak = sys.argv[1]
anamode = sys.argv[2]


"""
note: 
for repeating the fit for a given configuration, 
the remaining sources may be commented out below
"""

files.append("upsilonYieldWeighted_nominal.root")

#files.append("upsilonYieldWeighted_AccLo.root")
#files.append("upsilonYieldWeighted_AccHi.root")
#files.append("upsilonYieldWeighted_EtrkLo.root")
#files.append("upsilonYieldWeighted_EtrkHi.root")
##files.append("upsilonYieldWeighted_EmuidLo.root")
##files.append("upsilonYieldWeighted_EmuidHi.root")
##files.append("upsilonYieldWeighted_EtrigLo.root")
##files.append("upsilonYieldWeighted_EtrigHi.root")
#files.append("upsilonYieldWeighted_EtrecoHi.root")
#files.append("upsilonYieldWeighted_EtrecoLo.root")
#files.append("upsilonYieldWeighted_ptscaleLo.root")
#files.append("upsilonYieldWeighted_ptscaleHi.root")
#files.append("upsilonYieldWeighted_ptresoLo.root")
#files.append("upsilonYieldWeighted_ptresoHi.root")
#files.append("upsilonYieldWeighted_ptspec.root")
#files.append("upsilonYieldWeighted_vtxpos.root")
#files.append("upsilonYieldWeighted_nofsr.root")
#files.append("upsilonYieldWeighted_tnpmc.root")
#files.append("upsilonYieldWeighted_mctrue.root")
#files.append("upsilonYieldWeighted_tnpmcUps.root")
#files.append("upsilonYieldWeighted_linear812.root")
#files.append("upsilonYieldWeighted_otherLo.root")
#files.append("upsilonYieldWeighted_otherHi.root")

if peak == '1' :
    print "upsilon 1s\n"
    if pol:
        files.append("upsilonYieldWeighted_helT.root")
        files.append("upsilonYieldWeighted_helL.root")
        files.append("upsilonYieldWeighted_csT.root")
        files.append("upsilonYieldWeighted_csL.root")
        
elif peak == '2' :
    print "upsilon 2s\n"
    files.append("upsilonYieldWeighted_2s.root")
    if pol:
        files.append("upsilonYieldWeighted_csL_2s.root")
        files.append("upsilonYieldWeighted_csT_2s.root")
        files.append("upsilonYieldWeighted_helL_2s.root")
        files.append("upsilonYieldWeighted_helT_2s.root")
        
elif peak == '3' :
    print "upsilon 3s\n"
    files.append("upsilonYieldWeighted_3s.root")
    if pol:
        files.append("upsilonYieldWeighted_csL_3s.root")
        files.append("upsilonYieldWeighted_csT_3s.root")
        files.append("upsilonYieldWeighted_helL_3s.root")
        files.append("upsilonYieldWeighted_helT_3s.root")


print 'processing:'
print files


fit_cmd = "root -l -b -q oniafitter.C"

mode_dir = 'mode' + anamode
if not os.path.isdir(mode_dir):
    os.system('mkdir ' + mode_dir)

for fin in files : 
    ftmp = fin.replace('/',"_").replace(".root","")
    dres = mode_dir + '/' + "fitres_" + peak + "s_" + ftmp + '/'
    flog = dres + "log"
    fout = dres + "xsection.root"
    if fin.find("_linear")==-1 :
        bgflag = "false"
    else :
        bgflag = "true"
    command = fit_cmd + '\\(\\"' + fin + '\\",\\"' + fout + '\\",\\"' + dres + '\\",\\"' + peak + '\\",\\"' + bgflag + '\\",\\"' + 'mode' + anamode + '\\"\\) >& ' + flog + ' &'
    print command 
    #sys.exit(0)
    if not os.path.isdir(dres) :
        os.system('mkdir ' + dres)
    if os.path.isfile(flog) :
        os.system('mv ' + flog + ' ' + flog + ".old" )
    os.system(command)

    #time.sleep(120) #wait 2 min between fits 
    while True:
        time.sleep(180)        
        #if os.path.isfile(dres + 'mass_raw.pdf'):
        if os.path.isfile(dres + 'massfit_rap0_pt0.pdf'):
            break
        
        
