import platform
import os
import sys
import glob, shutil

def install_python_interface(project_root_dir, ocp_name, install_prefix=None):
    if install_prefix is None:
        python_version = 'python' + str(sys.version_info.major) + '.' + str(sys.version_info.minor)
        if platform.system() == 'Windows':
            install_prefix = os.path.join(os.path.abspath(os.environ['HOMEPATH']), '.local', 'lib', python_version, 'site-packages')
        else:
            install_prefix = os.path.join(os.path.abspath(os.environ['HOME']), '.local', 'lib', python_version, 'site-packages')
    install_destination = os.path.join(os.path.abspath(install_prefix), 'cgmres')
    build_dir = os.path.join(project_root_dir, 'build')
    pybind11_sharedlibs_dir = os.path.join(build_dir, 'python', ocp_name)
    pybind11_sharedlibs = glob.glob(os.path.join(pybind11_sharedlibs_dir, '*.so')) \
                            + glob.glob(os.path.join(pybind11_sharedlibs_dir, '*.dylib')) \
                            + glob.glob(os.path.join(pybind11_sharedlibs_dir, '*.pyd')) 
    pybind11_sharedlibs_common_dir = os.path.join(build_dir, 'python', 'common')
    pybind11_sharedlibs_common = glob.glob(os.path.join(pybind11_sharedlibs_common_dir, '*.so')) \
                                  + glob.glob(os.path.join(pybind11_sharedlibs_common_dir, '*.dylib')) \
                                  + glob.glob(os.path.join(pybind11_sharedlibs_common_dir, '*.pyd')) 
    print('Collected Python shared libs: ', pybind11_sharedlibs)
    print('Collected Python common shared libs: ', pybind11_sharedlibs_common)
    os.makedirs(os.path.join(install_destination, ocp_name), exist_ok=True)
    os.makedirs(os.path.join(install_destination, 'common'), exist_ok=True)
    for e in pybind11_sharedlibs:
        shutil.copy(e, os.path.join(install_destination, ocp_name))
    for e in pybind11_sharedlibs_common:
        shutil.copy(e, os.path.join(install_destination, 'common'))
    python_files = glob.glob(os.path.join(os.path.join(project_root_dir, 'python', ocp_name), '*.py'))
    python_files_common = glob.glob(os.path.join(os.path.join(project_root_dir, 'python', 'common'), '*.py'))
    for e in python_files:
        shutil.copy(e, os.path.join(install_destination, ocp_name))
    for e in python_files_common:
        shutil.copy(e, os.path.join(install_destination, 'common'))
    print('Collected Python files: ', python_files)
    print('Collected Python common files: ', python_files_common)
    print('\nPython interfaces have been installed at ' + str(install_prefix))
    print('To use Python interfaces, run \n')
    print('    export PYTHONPATH=$PYTHONPATH:' + str(install_prefix) + '\n')
    print('in the terminal to recognize the PYTHONPATH temporary.')
    print('Or set the PATH in Ubuntu as\n')
    print('    echo export PYTHONPATH=$PYTHONPATH:' + str(install_prefix) + ' >> ~/.bashrc\n')
    print('or in Mac OSX as\n')
    print('    echo export PYTHONPATH=$PYTHONPATH:' + str(install_prefix) + ' >> ~/.zshrc\n')


if __name__ == '__main__':
    assert len(sys.argv) >= 3
    project_root_dir = sys.argv[1]
    ocp_name = sys.argv[2]
    if len(sys.argv) >= 4:
        install_prefix = sys.argv[3]
    else:
        install_prefix = None
    install_python_interface(project_root_dir, ocp_name, install_prefix)