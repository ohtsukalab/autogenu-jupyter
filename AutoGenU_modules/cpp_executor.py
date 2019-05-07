import subprocess
import platform


def setCMake(simulation_name, MSYS=False):
    if(platform.system() == 'Windows'):
        subprocess.run(['mkdir', 'build'], cwd='models/'+simulation_name, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell=True)
        if(MSYS == True):
            proc = subprocess.Popen(['cmake', '../../..', '-G', 'MSYS Makefiles'], cwd='models/'+simulation_name+'/build', stdout = subprocess.PIPE, stderr = subprocess.STDOUT, shell=True)
        else:
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


def removeBuildDir(simulation_name):
    if(platform.system() == 'Windows'):
        subprocess.run(['rmdir', '/q', '/s', 'build'], cwd='models/'+simulation_name, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell=True)
    else:
        subprocess.run(['rm', '-r', 'build'], cwd='models/'+simulation_name, stdout = subprocess.PIPE, stderr = subprocess.PIPE)


def makeAndRun(simulation_name):
    if(platform.system() == 'Windows'):
        proc = subprocess.Popen(['cmake', '--build', '.'], cwd='models/'+simulation_name+'/build', stdout = subprocess.PIPE, stderr = subprocess.STDOUT, shell=True)
        for line in iter(proc.stdout.readline,b''):
            print(line.rstrip().decode("utf8"))
        print('\n')
        proc = subprocess.Popen(['main.exe'], cwd='models/'+simulation_name+'/build', stdout = subprocess.PIPE, stderr = subprocess.STDOUT, shell=True)
        for line in iter(proc.stdout.readline,b''):
            print(line.rstrip().decode("utf8"))
        subprocess.run(['rmdir', '/q', '/s', 'simulation_result'], cwd='models/'+simulation_name, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell=True)
        subprocess.run(['move', 'simulation_result', '../'], cwd='models/'+simulation_name+'/build', stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell=True)
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




