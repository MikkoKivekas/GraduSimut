CC = g++
LDFLAGS = `pkg-config --libs ibsimu-1.0.6dev`
CXXFLAGS = -Wall -g `pkg-config --cflags ibsimu-1.0.6dev`

simu: simu.o
	$(CC) -o simu simu.o $(LDFLAGS)

simu.o: simu.cpp
	$(CC) -c -o simu.o simu.cpp $(CXXFLAGS)

clean:
	$(RM) *~ *.o simu