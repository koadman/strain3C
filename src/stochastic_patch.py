#!/usr/bin/env python3
#
import sys

search_1 = "for (int i = 1; i <= L_hic; ++i)"
search_2 = "for (int i = 1; i <= L_local; ++i)"

for line in sys.stdin:
  print(line.rstrip())
  if line.find(search_1) >= 0:
    search_1 += "gobbledeegook" # has been used so make it unfindable
    print("                if(rand() > RAND_MAX / subsample) continue; // stochastic gradient descent")
  if line.find(search_2) >= 0:
    search_2 += "gobbledeegook" # has been used so make it unfindable
    print("                if(rand() > RAND_MAX / subsample) continue; // stochastic gradient descent")
