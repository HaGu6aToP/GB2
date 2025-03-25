
# filename = main
# filename = garray_search_test
# start:
# 	gcc $(filename).c -omain.out -I/usr/lib/x86_64-linux-gnu/glib-2.0/include -I/usr/include/glib-2.0 -I/usr/include/x86_64-linux-gnu/flint -L/usr/lib/x86_64-linux-gnu -lmpfr -lflint -lgmp -lglib-2.0

# all: basis_tools

# tools:
# 	gcc tools.c -c -I/usr/lib/x86_64-linux-gnu/glib-2.0/include -I/usr/include/glib-2.0 -I/usr/include/x86_64-linux-gnu/flint -L/usr/lib/x86_64-linux-gnu -lmpfr -lflint -lgmp -lglib-2.0

# basis_tools:
# 	gcc basis_tools -c -I/usr/lib/x86_64-linux-gnu/glib-2.0/include -I/usr/include/glib-2.0 -I/usr/include/x86_64-linux-gnu/flint -L/usr/lib/x86_64-linux-gnu -lmpfr -lflint -lgmp -lglib-2.0

# main: 
# 	gcc $(filename).c -c -I/usr/lib/x86_64-linux-gnu/glib-2.0/include -I/usr/include/glib-2.0 -I/usr/include/x86_64-linux-gnu/flint -L/usr/lib/x86_64-linux-gnu -lmpfr -lflint -lgmp -lglib-2.0

# start: tools basis_tools main
# 	gcc tools.o basis_tools.o main.o -omain

target = test
# src = $(wildcard *.c)
srcnames = tools.c basis_tools.c buchberger.c
srcdir = ../src/
incdir = /usr/code/GB2/include/
src = $(addprefix $(srcdir), $(srcnames)) min_find_test.c main.c
obj = $(patsubst %.c, %.o, $(src))



cflags = -I/usr/lib/x86_64-linux-gnu/glib-2.0/include \
		-I/usr/include/glib-2.0 \
		-I/usr/include/x86_64-linux-gnu/flint \
		-I/usr/code/GB2 \
		-I../include/

ldflags = -L/usr/lib/x86_64-linux-gnu \
		-lmpfr -lflint -lgmp -lglib-2.0

$(target) : $(obj)
	gcc $(obj) -o$(target) $(ldflags)

%.o : %.c 
	gcc -c $< -o $@ $(cflags) 

clean : 
	rm $(target) *.o