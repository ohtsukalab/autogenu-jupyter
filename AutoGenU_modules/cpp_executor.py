import subprocess
import platform


def setCMake(simulation_name):
    subprocess.run(['mkdir', 'build'], cwd='models/'+simulation_name, shell=True)
    subprocess.run(['cmake', '../../..'], cwd='models/'+simulation_name+'/build', shell=True)


def makeAndRun(simulation_name):
    if(platform.system() == 'Windows'):
        subprocess.run(['cmake', '--build', '.'], cwd='models/'+simulation_name+'/build', shell=True)
        subprocess.run(['a.exe'], cwd='models/'+simulation_name+'/build', shell=True)
    else:
        subprocess.run(['cmake', '--build', '.'], cwd='models/'+simulation_name+'/build', shell=True)
        subprocess.run(['./a.out'], cwd='models/'+simulation_name+'/build', shell=True)
    subprocess.run(['rm', '-r', 'simulation_result'], cwd='models/'+simulation_name, shell=True)
    subprocess.run(['mv', 'simulation_result', '../'], cwd='models/'+simulation_name+'/build', shell=True)




