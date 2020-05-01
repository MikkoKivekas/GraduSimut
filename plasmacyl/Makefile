CC = g++
LDFLAGS = `pkg-config --libs ibsimu-1.0.6dev`
CXXFLAGS = -Wall -g `pkg-config --cflags ibsimu-1.0.6dev`

plasmacyl: plasmacyl.o
	$(CC) -o plasmacyl plasmacyl.o $(LDFLAGS)

plasmacyl.o: plasmacyl.cpp
	$(CC) -c -o plasmacyl.o plasmacyl.cpp $(CXXFLAGS)

clean:
	$(RM) *~ *.o plasmacyl
