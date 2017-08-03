CC = g++
CFLAGS = -O3 -Wall -std=c++14 -Wno-deprecated-declarations -ftemplate-backtrace-limit=0 `pkg-config --cflags eigen3 opencv` -I./
LIBS = `pkg-config --libs eigen3 opencv`
OBJS = test.o
TARGET = test

$(TARGET): $(OBJS)
	$(CC) $^ -o $@ $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@
%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -rf *.o *.swp $(TARGET)
