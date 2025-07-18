import matplotlib.pyplot as plt
d = open("/Users/tomszwarcer/Documents/CERN/spacecharge-sim/build/ion_pos.txt", "r")
positions = d.readlines()
positions = [float(i) for i in positions]
plt.hist(positions,40)
plt.title("Ion positions")
plt.xlabel("y [cm]")
plt.ylabel("Counts")
plt.savefig("/Users/tomszwarcer/Documents/CERN/spacecharge-sim/build/ion_distribution.png")