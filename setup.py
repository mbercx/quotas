from setuptools import setup, find_packages

setup(
    name='quotas',
    version='0.1',
    packages=find_packages(exclude=['docs']),
    install_requires=[
        'scipy',
        'numpy',
        'click',
        'pymatgen'
    ],
    entry_points='''
        [console_scripts]
        quotas=quotas.cli.cli:main
    ''',
)