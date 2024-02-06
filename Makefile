CC = g++
CFLAGS = -O3 -g -DGL_GLEXT_PROTOTYPES -DGL_DO_NOT_WARN_IF_MULTI_GL_VERSION_HEADERS_INCLUDED -DOSX -Wno-deprecated-register -Wno-deprecated-declarations -Wno-shift-op-parentheses
INCFLAGS = -I./glm-0.9.7.1 -I/usr/X11/include -I./include/
LDFLAGS =  -L./lib/mac/ -L./lib/nix/ \
		-L"/usr/local/lib" \
		-lGL -lGLU -lm -lstdc++ -lfreeimage

RM = /bin/rm -f 
all: transforms
transforms: main.o Transform.o readfile.o variables.h readfile.h Transform.h
	$(CC) $(CFLAGS) -o transforms main.o Transform.o readfile.o  $(INCFLAGS) $(LDFLAGS)
main.o: main.cpp Transform.h variables.h
	$(CC) $(CFLAGS) $(INCFLAGS) -c main.cpp
readfile.o: readfile.cpp readfile.h variables.h 
	$(CC) $(CFLAGS) $(INCFLAGS) -c readfile.cpp
Transform.o: Transform.cpp Transform.h 
	$(CC) $(CFLAGS) $(INCFLAGS) -c Transform.cpp
clean: 
	$(RM) *.o transforms *.png


 
