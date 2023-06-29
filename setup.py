from setuptools import setup, find_packages
import TEster
setup(
    name='TEster',
    version=TEster.__version__,
    description='nester testing tool',
    author=TEster.__author__,
    packages=find_packages(),
    install_requires=[
        'bcbio-gff>=0.6.4',
        'biopython>=1.70',
        'click>=6.7',
        'networkx>=2.1',
        'PyYAML>=3.12',
        'ruamel.yaml>=0.16.10',
        'pyfastx>=0.6.6',
        'numpy>=1.17',
        'scipy>=1.4.1'
    ],
    entry_points={
        'console_scripts': [
            'TEster = TEster.toolchain.TEster:main'
        ]
    }
)

