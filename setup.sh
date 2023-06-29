#!/bin/sh


rm -rf dist build TEster.egg-info

if [ "$1" == "--user" ]
  then
	if  [ ! -f ~/.local/nested/config.yml ]
    	then
        	echo "Error: ~/.local/nested/config.yml file not found, try running TE-nester's setup.sh"
        	exit 1
    	fi
        python3 setup.py install --user

else
    if  [ ! -f /etc/nested/config.yml ]
    then
        echo "Error: /etc/nested/config.yml file not found, try running TE-nester's setup.sh"
        exit 1
    fi
    python3 setup.py install
fi
