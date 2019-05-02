import subprocess
import platform


def setCMake(simulation_name):
    if(platform.system() == 'Windows'):
        subprocess.run(['mkdir', 'build'], cwd='models/'+simulation_name, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell=True)
        proc = subprocess.Popen(['cmake', '../../..', '-G', 'MinGW Makefiles'], cwd='models/'+simulation_name+'/build', stdout = subprocess.PIPE, stderr = subprocess.STDOUT, shell=True)
        for line in iter(proc.stdout.readline,b''):
            print(line.rstrip().decode("utf8"))
        print('\n')
    else:
        subprocess.run(['mkdir', 'build'], cwd='models/'+simulation_name, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        proc = subprocess.Popen(['cmake', '../../..'], cwd='models/'+simulation_name+'/build', stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
        for line in iter(proc.stdout.readline,b''):
            print(line.rstrip().decode("utf8"))
        print('\n')


def makeAndRun(simulation_name):
    if(platform.system() == 'Windows'):
        proc = subprocess.Popen(['cmake', '--build', '.'], cwd='models/'+simulation_name+'/build', stdout = subprocess.PIPE, stderr = subprocess.STDOUT, shell=True)
        for line in iter(proc.stdout.readline,b''):
            print(line.rstrip().decode("utf8"))
        print('\n')
        proc = subprocess.run(['main.exe'], cwd='models/'+simulation_name+'/build', stdout = subprocess.PIPE, stderr = subprocess.STDOUT, shell=True)
        for line in iter(proc.stdout.readline,b''):
            print(line.rstrip().decode("utf8"))
        subprocess.run(['rd', '/Q', '/S', 'simulation_result'], cwd='models/'+simulation_name, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        subprocess.run(['move', 'simulation_result', '../'], cwd='models/'+simulation_name+'/build', stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    else:
        proc = subprocess.Popen(['cmake', '--build', '.'], cwd='models/'+simulation_name+'/build', stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
        for line in iter(proc.stdout.readline,b''):
            print(line.rstrip().decode("utf8"))
        print('\n')
        proc = subprocess.Popen(['./a.out'], cwd='models/'+simulation_name+'/build', stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
        for line in iter(proc.stdout.readline,b''):
            print(line.rstrip().decode("utf8"))
        subprocess.run(['rm', '-r', 'simulation_result'], cwd='models/'+simulation_name, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell=True)
        subprocess.run(['mv', 'simulation_result', '../'], cwd='models/'+simulation_name+'/build', stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell=True)




