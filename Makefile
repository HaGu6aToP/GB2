
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

target = main
# src = $(wildcard *.c)
src = main.c tools.c basis_tools.c buchberger.c
obj = $(patsubst %.c, %.o, $(src))



cflags = -I/usr/lib/x86_64-linux-gnu/glib-2.0/include \
		-I/usr/include/glib-2.0 \
		-I/usr/include/x86_64-linux-gnu/flint \
		-I/usr/code/GB2

ldflags = -L/usr/lib/x86_64-linux-gnu \
		-lmpfr -lflint -lgmp -lglib-2.0

$(target) : $(obj)
	gcc $(obj) -o$(target) $(ldflags)

%.o : %.c 
	gcc -c $< -o $@ $(cflags) 

clean : 
	rm $(target) *.o