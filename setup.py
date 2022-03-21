from setuptools import setup

# get the version here
pkg_vars = {}

with open("version.py") as fp:
    exec(fp.read(), pkg_vars)

setup(
    name='ztf_pipeutils',
    version=pkg_vars['__version__'],
    description='',
    url='http://github.com/pgris/ztf_pipeutils',
    author='M. Bailleul, Ph. Gris',
    author_email='manon.bailleul@etu.uca.fr,philippe.gris@clermont.in2p3.fr',
    license='BSD',
    packages=['ztf_pipeutils'],
    python_requires='>=3.5',
    zip_safe=False,
)
