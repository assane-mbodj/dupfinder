#!/bin/bash
set -e

#lumpy-sv
sudo apt install lumpy-sv && \

#smoove
wget https://github.com/brentp/smoove/releases/download/v0.2.8/smoove && \
chmod 755 smoove
mv smoove ~/anaconda3/envs/dupfinder_env/bin/ && \

#Duphold
wget https://github.com/brentp/duphold/releases/download/v0.2.3/duphold 

chmod 755 duphold && \
mv duphold ~/anaconda3/envs/dupfinder_env/bin/ 

rm -rf lumpy-sv
