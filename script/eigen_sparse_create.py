# string = '''1 0 0 0 0 0 0 1 1 4 0 0
# 1 6 0 2 0 0 0 0 0 0 0 0
# 0 1 6 1 1 1 4 0 0 0 0 0
# 0 0 0 0 0 1 0 0 0 0 6 2'''

string = '''1 0 1 0 0
1 1 0 0 0
0 1 0 1 0
0 0 0 6 1'''

lines = string.split('\n')

for (i, line) in enumerate(lines):
    l = line.strip()
    j = 0
    for char in l:
        if char == ' ':
            continue
        if char != '0':
            print(f"M.insert({i}, {j}) = {char};")
        j += 1