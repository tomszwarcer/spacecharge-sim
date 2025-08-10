out = open("all_runs.txt", "w")
dv = [550]
counter = -1
num_runs = 400
start = 0
for i in range(start,start + num_runs):
    if i%(num_runs/len(dv)) == 0:
        counter += 1
    out.write(str(i)+","+str(dv[counter]) + "\n")
out.close