
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


# target = main
# # src = $(wildcard *.c)
# srcnames = main.c tools.c basis_tools.c buchberger.c f4.c
# srcdir = ./src/
# incdir = ./include/
# src = $(addprefix $(srcdir), $(srcnames))
# obj = $(patsubst %.c, %.o, $(src))

# cflags = -I/usr/lib/x86_64-linux-gnu/glib-2.0/include \
# 		-I/usr/include/glib-2.0 \
# 		-I/usr/include/x86_64-linux-gnu/flint \
# 		-I/usr/code/GB2 \
# 		-I./include/ \

# ldflags = -L/usr/lib/x86_64-linux-gnu \
# 		-lmpfr -lflint -lgmp -lglib-2.0

# $(target) : $(obj)
# 	g++ $(obj) -o$(target) $(ldflags)

# %.o : %.c 
# 	g++ -xc -c $< -o $@ $(cflags) 

# clean : 
# 	rm $(target) *.o


# Определение цели
target = main

# Получение списка исходных файлов C и C++
src_c =  $(wildcard ./src/*.c)
src_cpp = $(wildcard ./src/*.cpp)
srcnames = $(notdir $(src_c)) $(notdir $(src_cpp))
srcdir = ./src/
incdir = ./include/
src = $(addprefix $(srcdir), $(srcnames))

# Получение списка объектных файлов
obj_c = $(patsubst %.c, %.o, $(src_c))
obj_cpp = $(patsubst %.cpp, %.o, $(src_cpp))
obj = $(obj_c) $(obj_cpp)

# Флаги компиляции
cflags = -I/usr/lib/x86_64-linux-gnu/glib-2.0/include \
		-I/usr/include/glib-2.0 \
		-I/usr/include/x86_64-linux-gnu/flint \
		-I/usr/code/GB2 \
		-I../SparseRREFF \
		-I$(incdir)

# Флаги линковки
ldflags = -L/usr/lib/x86_64-linux-gnu \
		-lmpfr -lflint -lgmp -lglib-2.0

# Правила компиляции
$(target) : $(obj)
	g++ $(obj) -o $(target) $(ldflags)

%.o : %.cpp
	g++ -fpermissive -c $< -o $@ $(cflags)

%.o : %.c
	g++ -fpermissive -c $< -o $@ $(cflags)

# Правило для очистки
clean:
	rm -f $(target) $(obj)