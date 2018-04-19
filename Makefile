CXXFLAGS = -O2 -Wall -Wno-deprecated -DNUMBITS=32

TF_CFLAGS=$(shell python -c 'import tensorflow as tf; print(" ".join(tf.sysconfig.get_compile_flags()))')
TF_LFLAGS=$(shell python -c 'import tensorflow as tf; print(" ".join(tf.sysconfig.get_link_flags()))')
INCLUDES =  -I./multilevelLib -I./metisLib -I/usr/local/lib -I/opt/local/include -L./multilevelLib -L./metisLib

LIBS = -lmultilevel -lmetis -lm

default:
	(cd metisLib ; make )
	(cd multilevelLib ; make )
	g++ -std=c++11 -shared graclus.cc -o graclus.so -fPIC ${INCLUDES} ${TF_CFLAGS} ${TF_LFLAGS} -O2

clean:
	(cd metisLib ; make clean )
	(cd multilevelLib ; make clean )
	(cd programs ; make clean )

realclean:
	(cd metisLib ; make realclean )
	(cd multilevelLib ; make realclean )
	(cd programs ; make realclean )
	(rm *.a)