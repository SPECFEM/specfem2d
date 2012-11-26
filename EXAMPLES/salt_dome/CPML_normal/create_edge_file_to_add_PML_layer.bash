#!/bin/bash

gawk '{print $3  "\n" $4}' a | sort -n | uniq > bottom_edge

