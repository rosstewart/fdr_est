#!/bin/bash

nohup matlab -nodisplay -r "$1" > matlab$2.log 2>&1 &


