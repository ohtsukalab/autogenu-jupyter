import subprocess
import platform


def setCMake(simulation_name):
    subprocess.run(['mkdir', 'build'], cwd='models/'+simulation_name)
    subprocess.run(['cmake', '../../..'], cwd='models/'+simulation_name+'/build')


def makeAndRun(simulation_name):
    if(platform.system() == 'Windows'):
        subprocess.run(['cmake', '--build', '.'], cwd='models/'+simulation_name+'/build')
        subprocess.run(['a.exe'], cwd='models/'+simulation_name+'/build')
    else:
        subprocess.run(['cmake', '--build', '.'], cwd='models/'+simulation_name+'/build')
        subprocess.run(['./a.out'], cwd='models/'+simulation_name+'/build')
    subprocess.run(['rm', '-r', 'simulation_result'], cwd='models/'+simulation_name)
    subprocess.run(['mv', 'simulation_result', '../'], cwd='models/'+simulation_name+'/build')




