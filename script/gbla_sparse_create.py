# string = '''1 0 0 0 0 0 0 1 1 4 0 0
# 1 6 0 2 0 0 0 0 0 0 0 0
# 0 1 6 1 1 1 4 0 0 0 0 0
# 0 0 0 0 0 1 0 0 0 0 6 2'''

string = ""

f = open("aboba.txt", "r")
string = f.read()

# string = '''0 1 0 0 1 0 0 0 0 0 0 
# 0 1 6 0 0 0 0 0 0 0 0 
# 0 0 0 0 0 6 0 0 0 1 0 
# 0 0 0 0 0 1 0 6 0 0 0 
# 1 0 0 1 0 0 0 0 0 0 0 
# 6 0 1 0 0 0 0 0 0 0 0 
# 0 0 1 0 0 0 1 0 0 0 0 
# 0 0 0 6 0 0 1 0 0 0 0 
# 0 0 0 0 6 0 0 0 1 0 0 
# 0 0 0 0 0 0 6 0 0 0 1 
# 0 0 0 0 0 0 0 0 1 0 6'''


# string = '''1 0 1 0 0
# 1 1 0 0 0
# 0 1 0 1 0
# 0 0 0 6 1'''

lines = string.split('\n')
if lines[-1] == '':
    lines.pop()
pos = ''
f = open("res.txt", "w")
mod = 7
m = len(lines)
n = 0
for el in lines[0]:
    if el != ' ':
        n+=1


f.write(f"{mod} {m} {n}\n")

for line in lines:
    count = 0
    j = 0
    data = ''
    pos = ''
    for el in line:
        if el == ' ':
            continue
        if el == '0':
            j+=1
            continue
        pos += str(j) + ' '
        data += str(el) + ' '
        count += 1
        j += 1

    f.write(str(count) + '\n')
    f.write(pos + '\n')
    f.write(data + '\n')


f.close()
