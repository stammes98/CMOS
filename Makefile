LFLAGS = -lueye_api -lcfitsio

make: main

main: main.o CMOSUtils.o
	g++ main.o CMOSUtils.o -o main $(LFLAGS)

main.o: main.cpp
	g++ -c main.cpp -o main.o -w
CMOSUtils.o: CMOSUtils.cpp
	g++ -c CMOSUtils.cpp -o CMOSUtils.o -w

clean:
	rm -r *.o

run:
	ulimit -s unlimited
	./main
