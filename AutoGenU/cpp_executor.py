import subprocess
import platform


def set_cmake(model_name, MSYS=False):
    """ Makes a directory for build and generate Makefiles using CMake.

        Args: 
            model_name: A string representing the name of the simulation model.
            MSYS: An optional variable, which should be true when you use 
                Windows and MSYS is default compiler.
    """
    if(platform.system() == 'Windows'):
        subprocess.run(
            ['mkdir', 'build'], 
            cwd='models/'+model_name, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, 
            shell=True
        )
        if(MSYS == True):
            proc = subprocess.Popen(
                ['cmake', '../../..', '-G', 'MSYS Makefiles'], 
                cwd='models/'+model_name+'/build', 
                stdout=subprocess.PIPE, 
                stderr=subprocess.STDOUT, 
                shell=True
            )
        else:
            proc = subprocess.Popen(
                ['cmake', '../../..', '-G', 'MinGW Makefiles'], 
                cwd='models/'+model_name+'/build', 
                stdout=subprocess.PIPE, 
                stderr=subprocess.STDOUT, 
                shell=True
            )
        for line in iter(proc.stdout.readline, b''):
            print(line.rstrip().decode("utf8"))
        print('\n')
    else:
        subprocess.run(
            ['mkdir', 'build'], 
            cwd='models/'+model_name, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE
        )
        proc = subprocess.Popen(
            ['cmake', '../../..'], 
            cwd='models/'+model_name+'/build', 
            stdout=subprocess.PIPE, 
            stderr=subprocess.STDOUT
        )
        for line in iter(proc.stdout.readline, b''):
            print(line.rstrip().decode("utf8"))
        print('\n')


def remove_build_dir(model_name):
    """ Removes a build directory. This function is mainly for Windows users 
        with MSYS.

        Args: 
            model_name: A string representing the name of the simulation model.
    """
    if(platform.system() == 'Windows'):
        subprocess.run(
            ['rmdir', '/q', '/s', 'build'], 
            cwd='models/'+model_name, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, 
            shell=True
        )
    else:
        subprocess.run(
            ['rm', '-r', 'build'],
            cwd='models/'+model_name, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE
        )


def make_and_run(model_name):
    """ Generates an execute file of a simulation and run it.

        Args: 
            model_name: A string representing the name of the simulation model.
    """
    if(platform.system() == 'Windows'):
        proc = subprocess.Popen(
            ['cmake', '--build', '.'], 
            cwd='models/'+model_name+'/build', 
            stdout=subprocess.PIPE, 
            stderr=subprocess.STDOUT, 
            shell=True
        )
        for line in iter(proc.stdout.readline,b''):
            print(line.rstrip().decode("utf8"))
        print('\n')
        subprocess.run(
            ['rmdir', '/q', '/s', 'simulation_result'], 
            cwd='models/'+model_name, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, 
            shell=True
        )
        proc = subprocess.Popen(
            ['main.exe'], 
            cwd='models/'+model_name+'/build', 
            stdout=subprocess.PIPE, 
            stderr=subprocess.STDOUT, 
            shell=True
        )
        for line in iter(proc.stdout.readline, b''):
            print(line.rstrip().decode("utf8"))
    else:
        proc = subprocess.Popen(
            ['cmake', '--build', '.'], 
            cwd='models/'+model_name+'/build', 
            stdout = subprocess.PIPE, 
            stderr = subprocess.STDOUT
        )
        for line in iter(proc.stdout.readline, b''):
            print(line.rstrip().decode("utf8"))
        print('\n')
        subprocess.run(
            ['rm', '-rf', 'simulation_result'], 
            cwd='models/'+model_name, 
            stdout = subprocess.PIPE, 
            stderr = subprocess.PIPE, 
            shell=True
        )
        proc = subprocess.Popen(
            ['./a.out'], 
            cwd='models/'+model_name+'/build', 
            stdout=subprocess.PIPE, 
            stderr=subprocess.STDOUT
        )
        for line in iter(proc.stdout.readline, b''):
            print(line.rstrip().decode("utf8"))