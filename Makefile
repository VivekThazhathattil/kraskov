CC = h5pcc #h5cc
MATH_FLAG = -lm
OPT_FLAG = -O3
PARALLEL_FLAG = -fopenmp
SRC=src
OUT_DIR=out
OUT_NAME=kraskov
INCLUDE=include/

kraskov:
	mkdir -p $(OUT_DIR);
	$(CC) $(MATH_FLAG) \
		$(OPT_FLAG) $(PARALLEL_FLAG)\
		$(SRC)/oth_utils.c \
		$(SRC)/mat_utils.c \
		$(SRC)/digamma.c \
		$(SRC)/data_handler.c \
		$(SRC)/print_utils.c \
		$(SRC)/kraskov.c \
		$(SRC)/main.c \
		-o $(OUT_NAME) -I$(INCLUDE);
	mv *.o $(OUT_DIR);

clean:
	rm -rf $(OUT_DIR);
	rm $(OUT_NAME);
