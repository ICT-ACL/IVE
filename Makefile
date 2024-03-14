# compile main.cpp to main in the current directory

# The compiler: g++ for C++ program, gcc for C program

CC = g++

# The build target executable:

TARGET = main
FLAGS =  -std=c++11 -O3 -Wall -march=native

all: $(TARGET)

$(TARGET): main.cpp
	$(CC) $(FLAGS) -o $(TARGET) main.cpp

clean:
	$(RM) $(TARGET)
