#!/bin/bash
#gcc -lm digamma.c main.c
h5pcc -lm mat_utils.c digamma.c data_handler.c print_utils.c main.c
