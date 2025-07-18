import imageio.v2 as iio
import os

field = True


d = open("/Users/tomszwarcer/Documents/CERN/spacecharge-sim/build/filenames_drift.txt", "r")

filenames_d = d.readlines()

s = "/Users/tomszwarcer/Documents/CERN/spacecharge-sim/build/"
for i in range(len(filenames_d)):
    filenames_d[i] = s + filenames_d[i]
    filenames_d[i] = filenames_d[i][:-1]

w = iio.get_writer("/Users/tomszwarcer/Documents/CERN/spacecharge-sim/movies/drift_dynamic.gif", format="gif", mode='I', fps=10, loop=100)
for frame in filenames_d:
    w.append_data(iio.imread(frame))
w.close()

d.close()

if field:
    f = open("/Users/tomszwarcer/Documents/CERN/spacecharge-sim/build/filenames_field.txt", "r")
    filenames_f = f.readlines()
    for i in range(len(filenames_f)):
        filenames_f[i] = s + filenames_f[i]
        filenames_f[i] = filenames_f[i][:-1]
    x = iio.get_writer("/Users/tomszwarcer/Documents/CERN/spacecharge-sim/movies/field_dynamic.gif", format="gif", mode='I', fps=10, loop=100)
    for frame in filenames_f:
        x.append_data(iio.imread(frame)) 
    x.close()

    f.close()