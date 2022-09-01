import platform
import os
import sys
import glob, shutil

def install_python_interface(ocp_name, install_prefix=None):
    if install_prefix is None:
        python_version = 'python' + str(sys.version_info.major) + '.' + str(sys.version_info.minor)
        if platform.system() == 'Windows':
            install_prefix = os.path.join(os.path.abspath(os.environ['HOMEPATH']), '.local/lib', python_version, 'site-packages')
        else:
            install_prefix = os.path.join(os.path.abspath(os.environ['HOME']), '.local/lib', python_version, 'site-packages')
        install_destination = os.path.join(os.path.abspath(install_prefix), 'cgmres')
    else:
        install_destination = os.path.join(os.path.abspath(install_prefix), 'cgmres')
    pybind11_sharedlibs = glob.glob('build/python/'+ocp_name+'/*.so') \
                            + glob.glob('build/python/'+ocp_name+'/*.dylib') \
                            + glob.glob('build/python/'+ocp_name+'/*.pyd')
    pybind11_sharedlibs_common = glob.glob('build/python/common/*.so') \
                                    + glob.glob('build/python/common/*.dylib') \
                                    + glob.glob('build/python/common/*.pyd')
    os.makedirs(os.path.join(install_destination, ocp_name), exist_ok=True)
    os.makedirs(os.path.join(install_destination, 'common'), exist_ok=True)
    for e in pybind11_sharedlibs:
        shutil.copy(e, str(os.path.join(install_destination, ocp_name)))
    for e in pybind11_sharedlibs_common:
        shutil.copy(e, str(os.path.join(install_destination, 'common')))
    python_files = glob.glob('python/'+ocp_name+'/*.py')
    python_files_common = glob.glob('python/common/*.py')
    for e in python_files:
        shutil.copy(e, str(os.path.join(install_destination, ocp_name)))
    for e in python_files_common:
        shutil.copy(e, str(os.path.join(install_destination, 'common')))
    print('Python interfaces have been installed at ' + str(install_prefix))
    print('To use Python interfaces, run \n')
    print('    export PYTHONPATH=$PYTHONPATH:' + str(install_prefix) + '\n')
    print('in the terminal to recognize the PYTHONPATH temporary.')
    print('Or set the PATH in Ubuntu as\n')
    print('    echo export PYTHONPATH=$PYTHONPATH:' + str(install_prefix) + ' >> ~./bashrc\n')
    print('or in Mac OSX as\n')
    print('    echo export PYTHONPATH=$PYTHONPATH:' + str(install_prefix) + ' >> ~./zshrc\n')


if __name__ == '__main__':
    assert len(sys.argv) >= 2
    ocp_name = sys.argv[1]
    install_prefix = None
    if len(sys.argv) >= 3:
        install_prefix = sys.argv[2]
    install_python_interface(ocp_name, install_prefix)
    