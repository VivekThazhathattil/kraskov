#!/bin/bash
#gcc -lm digamma.c main.c
#h5cc -lm mat_utils.c digamma.c data_handler.c print_utils.c main.c -o kraskov
h5pcc -lm src/mat_utils.c src/digamma.c src/data_handler.c src/print_utils.c src/main.c -o kraskov -Iinclude/

mv *.o out/
