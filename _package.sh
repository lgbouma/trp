#!/bin/bash

# define the py38 venv from scratch; use the latest trp version on disk
cd ../ ;
tar czf trp-0.0.0.tar.gz trp ;
mv trp-0.0.0.tar.gz $HOME/environments/ ;
cd $HOME/environments/ ;

# NOTE: one dependency is not included in this auto-updated: the
# complexrotators package, https://github.com/lgbouma/cpv, because I have not
# put it on pip yet.
pip3 install --target=$HOME/environments/py38 -r requirements.txt ;
pip3 install --target=$HOME/environments/py38 --upgrade trp-0.0.0.tar.gz ;
rm trp-0.0.0.tar.gz ;

# move tarballed environment to runtime directory for exporting on OSG
cd $HOME/environments;
tar czf py38.tar.gz py38

cd $HOME/proj/trp/drivers
mv $HOME/environments/py38.tar.gz .
