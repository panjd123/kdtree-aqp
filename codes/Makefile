all: libaqp.so

libaqp.so: libaqp.cc libaqp.h
	g++ -Ofast -shared -fPIC -o libaqp.so libaqp.cc

clean:
	rm -f libaqp.so