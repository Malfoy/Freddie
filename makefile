CC=g++
CFLAGS= -Wall -Ofast -std=c++11  -flto -pipe -funit-at-a-time -fopenmp -lz -Isparsepp -flto
LDFLAGS=-flto -lpthread -fopenmp -lz  -Isparsepp  -flto
EXEC=Freddie Jason


all: $(EXEC)


Freddie:   freddie.o
	$(CC) -o $@ $^ $(LDFLAGS)

freddie.o: freddie.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

Jason: Jason.o
	$(CC) -o $@ $^ $(LDFLAGS)

Jason.o: jason.cpp
	$(CC) -o $@ -c $< $(CFLAGS)


clean:
	rm -rf *.o
	rm -rf $(EXEC)


rebuild: clean $(EXEC)
