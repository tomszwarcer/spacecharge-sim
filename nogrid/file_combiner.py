from os import listdir
from os.path import isfile, join
onpath = "/afs/cern.ch/work/t/tszwarce/nogrid/gain_output/on/"
onfiles = ["/afs/cern.ch/work/t/tszwarce/nogrid/gain_output/on/" + f for f in listdir(onpath) if isfile(join(onpath, f))]
offpath = "/afs/cern.ch/work/t/tszwarce/nogrid/gain_output/off/"
offfiles = ["/afs/cern.ch/work/t/tszwarce/nogrid/gain_output/off/" +f for f in listdir(offpath) if isfile(join(offpath, f))]
outfileoff = open("/afs/cern.ch/work/t/tszwarce/nogrid/gain_output/size_no_sc_combined.txt","w")
outfileon = open("/afs/cern.ch/work/t/tszwarce/nogrid/gain_output/size_sc_combined.txt","w")

for i in onfiles:
    file = open(i)
    l = file.readlines()
    try:
        outfileon.write(l[0])
    except IndexError:
        pass
    file.close()

for i in offfiles:
    file = open(i)
    l = file.readlines()
    try:
        outfileoff.write(l[0])
    except IndexError:
        pass
    file.close()