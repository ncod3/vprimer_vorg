"""Minimal setup file for vprimer project."""

from setuptools import setup, find_packages

setup(
    name='vprimer',
    version='0.0.1',
    license='GPL',
    description='Vprimer',

    author='satoshi-natsume',
    author_email="s-natsume@ibrc.or.jp",
    url="https://github.com/ncod3/vprimer",

    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    entry_points={
        'console_scripts': ['vprimer = src.main:main']
        }
)

