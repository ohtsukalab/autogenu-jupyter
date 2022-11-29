from setuptools import setup


def _requires_from_file(filename):
    return open(filename).read().splitlines()


setup(
    name="autogenu-jupyter",
    version="1.0.0",
    description="An automatic code generator for nonlinear model predictive control (NMPC) and the continuation/GMRES method (C/GMRES) based numerical solvers for NMPC",
    license="MIT",
    author="Sotaro Katayama",
    url="https://github.com/ohtsukalab/autogenu-jupyter",
    packages=["autogenu"],
    include_package_data=True,
    zip_safe=False,
    install_requires=_requires_from_file('requirements.txt'),
)