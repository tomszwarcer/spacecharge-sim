import imageio.v2 as iio
import os

field = True


from os import listdir
from os.path import isfile, join
mypath = "/Users/tomszwarcer/Documents/CERN/spacecharge-sim/nogrid/frames/"
filenames_d = [f for f in listdir(mypath) if isfile(join(mypath, f))]

ffiles = [i for i in filenames_d if i[0]=='f']
dfiles = [i for i in filenames_d if i[0]=='a']

nums = [int(i[-7:-4]) for i in filenames_d if i[0]=='a']
dictionary = {}
for i in range(len(ffiles)):
    txt = int(ffiles[i][-7:-4])
    for j in range(len(dfiles)):
        if int(dfiles[j][-7:-4]) == txt:
            dictionary[txt] = [ffiles[i],dfiles[j]]

d = dict(sorted(dictionary.items()))

w = iio.get_writer("/Users/tomszwarcer/Documents/CERN/spacecharge-sim/movies/drift_dynamic.gif", format="gif", mode='I', fps=10, loop=100)
for i in d:
    w.append_data(iio.imread(mypath+str(d[i][1])))
w.close()

w = iio.get_writer("/Users/tomszwarcer/Documents/CERN/spacecharge-sim/movies/field_dynamic.gif", format="gif", mode='I', fps=10, loop=100)
for i in d:
    w.append_data(iio.imread(mypath+str(d[i][0])))
w.close()
