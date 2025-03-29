# Compiler and flags
CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++17

# Source files
SRCS = main.cpp Nonlinear_sol_Newton.cpp AnalyticSolver.cpp loadTofile.cpp Tube.cpp

# Object files (automatically generated from source files)
OBJS = $(SRCS:.cpp=.o)

# Executable name
TARGET = program

# Default rule: build the executable
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET)

# Rule to compile source files into object files
%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: clean