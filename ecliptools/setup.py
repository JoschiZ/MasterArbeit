#!/usr/bin/python
from setuptools import setup

setup(
    name='eclup-tools',
    version='0.1.0',
    py_modules=['yourscript'],
    install_requires=[
        'Click',
    ],
    entry_points={
        'console_scripts': [
            'test = yourscript:cli',
        ],
    },
)
