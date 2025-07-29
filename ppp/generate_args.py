out = open("all_runs.txt", "w")
dv = [480,500,520,540,560]
counter = -1
for i in range(100):
    if i%20 == 0:
        counter += 1
    out.write(str(i)+","+str(dv[counter]) + "\n")
out.close