CC = h5pcc #h5cc
FLAGS = -lm
SRC=src
OUT_DIR=out
OUT_NAME=kraskov
INCLUDE=include/

kraskov:
	mkdir -p $(OUT_DIR);
	$(CC) $(FLAGS) \
		$(SRC)/oth_utils.c \
		$(SRC)/mat_utils.c \
		$(SRC)/digamma.c \
		$(SRC)/data_handler.c \
		$(SRC)/print_utils.c \
		$(SRC)/kraskov.c \
		$(SRC)/main.c \
		-o $(OUT_NAME) -I$(INCLUDE);
	mv *.o $(OUT_DIR);
