CC = g++
CFLAGS = -g -Wall -std=c++11
TARGET = iubfc

all: $(TARGET)

$(TARGET) : $(TARGET).cpp
	$(CC) $(CFLAGS) -o $(TARGET) $(TARGET).cpp
