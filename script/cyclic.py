from sys import argv

_, n = argv
n = int(n)
n += 1

res = ""
res += str(n) + "\n" + str(7) + "\n" + str(n) + "\n"
for i in range(n-1):
    res += "x" + str(i) + " "
res += "x" + str(n-1)
res += "\n"

for i in range(n-1):
    res += "x" + str(i) + "+"

res += "x" + str(n-1) + "\n"

for i in range(n-1):
    res += "x" + str(i) + "*x" + str(i+1) + "+"

res += "x" + str(n-1) + "*x" + str(0) + "\n"

for i in range(3, n):
    for j in range(n-1):
        for k in range(i-1):
            res += "x" + str((j+k)%n) + "*"
        res += "x" + str((j+i-1)%n)
        res += "+"
    
    for k in range(i-1):
        res += "x" + str((n-1+k)%n) + "*"
    res += "x" + str((n-1 + i - 1)%n)
    res += "\n"

for i in range(n-1):
    res += "x" + str(i) + "*"
res += "x" + str(n-1)

res += "-1"

print(res)
