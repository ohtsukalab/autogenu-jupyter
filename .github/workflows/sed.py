files = ['cartpole.py', 'hexacopter.py', 'mobilerobot.py', 'pendubot.py']

for file in files:
    with open(file, "r") as f:
        lines = f.readlines()
    with open(file, "w") as f:
        for line in lines:
            if not "get_ipython()" in line:
                f.write(line)
