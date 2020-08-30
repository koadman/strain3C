#!/usr/bin/env python3
import json
import sys

json_file_in = open(sys.argv[1], 'r')
stan_json = json.load(json_file_in)
stan_json['K']=int(sys.argv[2])
json_file_out = open(sys.argv[3], 'w')
json.dump(stan_json, json_file_out)
