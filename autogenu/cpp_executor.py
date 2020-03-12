import subprocess
import platform


def set_cmake(model_name, generator='Auto', remove_build_dir=False):
    """ Makes a directory for build and generate Makefiles using CMake.

        Args: 
            model_name: A string representing the name of the simulation model.
            generator: An optional variable for Windows user to choose the
                generator. If 'MSYS', then 'MSYS Makefiles' is used. If 'MinGW',
                then 'MinGW Makefiles' is used. The default value is 'Auto' and
                the generator is selected automatically. If sh.exe exists in 
                your PATH, MSYS is choosed, and otherwise MinGW is used. If 
                different value from 'MSYS' and 'MinGW', generator is selected
                automatically.
            remove_build_dir: If true, the existing build directory is removed 
                and if False, the build directory is not removed.
                Need to be set True is you change CMake configuration, e.g., if 
                you change the generator. The default value is False.
    """
    if platform.system() == 'Windows':
        subprocess.run(
            ['mkdir', 'build'], 
            cwd='models/'+model_name, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, 
            shell=True
        )
        if generator == 'MSYS':
            proc = subprocess.Popen(
                ['cmake', '../../..', '-G', 'MSYS Makefiles'], 
                cwd='models/'+model_name+'/build', 
                stdout=subprocess.PIPE, 
                stderr=subprocess.STDOUT, 
                shell=True
            )
            for line in iter(proc.stdout.readline, b''):
                print(line.rstrip().decode("utf8"))
            print('\n')
        elif generator == 'MinGW':
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
            proc = subprocess.Popen(
                ['where', 'sh.exe'], 
                cwd='C:', 
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE, 
                shell=True
            )
            if proc.stderr.readline() == b'':
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
    if platform.system() == 'Windows':
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
    if platform.system() == 'Windows':
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