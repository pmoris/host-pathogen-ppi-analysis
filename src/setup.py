from setuptools import setup

setup(
    name='phppipy',
    version='1.0',
    description='Tools for analysing protein-protein interactions',
    url='https://github.com/pmoris/host-pathogen-ppi-analysis',
    author='Pieter Moris',
    author_email='pieter.moris@uantwerpen.be',
    license='MIT',
    packages=['phppipy'],
    install_requires=[
        'pandas', 'numpy', 'goscripts', 'statsmodels', 'scipy'
    ],
    zip_safe=False)
