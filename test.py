f = open("build/sorted_hap_ibd")
h = f.readlines()
f.close()
f = open("build/sorted_output")
m = f.readlines()
f.close()
m = list(set(m))
mset = set(m)
hset = set(h)
print(mset / hset)
print(len(mset & hset))

m = [i.strip().split() for i in m]
h = [i.strip().split() for i in h]


m = sorted(m, key = lambda x: float(x[7]), reverse=True)
# print(len(h))
# print(len(m))