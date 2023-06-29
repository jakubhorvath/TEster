# TEster

This package requires the base program TE-Greedy-Nester with which it is included to be installed first

## Dependencies 

    nested-nester
    nested-generator

## Python module dependencies 
    bcbio-gff >= 0.6.4
    biopython >= 1.70
    click >= 6.7
    networkx >= 2.1
    PyYAML >= 3.12
    ruamel.yaml >= 0.16.10
    pyfastx >= 0.6.6
    numpy >= 1.17
    scipy >= 1.4.1

## Installation

    $ chmod +x setup.sh
    
    # for the root privileged user 
    $ sudo ./setup.sh
    
    # or alternatively for user-space installation
    $ ./setup.sh --user


Keep in mind that when installing as root, TE-Greedy-Nester needs to be installed with root privileges as well and its <strong>config.yml</strong> file needs to be located in the <strong>/etc/nested/</strong> directory 

When installing as a regular user, the <strong>config.yml</strong> file needs to be located in the <strong>~/.local/etc/nested</strong> directory 

