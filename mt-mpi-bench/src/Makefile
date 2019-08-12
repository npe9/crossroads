CC=cc
CFLAGS=-O3 -std=c99
THREADED_CFLAGS=-DMTBENCH_THREADED -fopenmp -craympich-mt

all: mt-bench-pingpong \
	 st-bench-pingpong \
	 mt-bench-one-sided-passive \
	 st-bench-one-sided-passive \
	 mt-bench-msgrate \
     st-bench-msgrate

mt-bench-one-sided-passive: mt-onesided-passive.c
	$(CC) $(CFLAGS) $(THREADED_CFLAGS) -o mt-bench-one-sided-passive mt-onesided-passive.c

st-bench-one-sided-passive: mt-onesided-passive.c
	$(CC) $(CFLAGS) -o st-bench-one-sided-passive mt-onesided-passive.c

mt-bench-pingpong: mt-pingpong.c
	$(CC) $(CFLAGS) $(THREADED_CFLAGS) -o mt-bench-pingpong mt-pingpong.c

st-bench-pingpong: mt-pingpong.c
	$(CC) $(CFLAGS) -o st-bench-pingpong mt-pingpong.c

mt-bench-msgrate: mt-msgrate.c
	$(CC) $(CFLAGS) $(THREADED_CFLAGS) -o mt-bench-msgrate mt-msgrate.c

st-bench-msgrate: mt-msgrate.c
	$(CC) $(CFLAGS) -o st-bench-msgrate mt-msgrate.c
	

clean:
	rm mt-bench-pingpong \
		st-bench-pingpong \
		mt-bench-one-sided-passive \
		st-bench-one-sided-passive \
		mt-bench-msgrate \
		st-bench-msgrate
