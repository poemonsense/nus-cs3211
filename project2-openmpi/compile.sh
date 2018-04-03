mpicc -o pool pool.c pfile.c spec.c mympi.c logging.c --std=c99 -Wall -Werror -lpthread -lm $*
