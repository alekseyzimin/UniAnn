# Compiler and flags
CXX      := g++
CXXFLAGS := -std=c++17 -O3 -Wall -Wextra -march=native
LDFLAGS  := 

# Target
TARGET := uniann 
SRC    := uniann.cpp
OBJ    := $(SRC:.cpp=.o)

# Default rule
all: $(TARGET)

# Link
$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# Compile with automatic dependency generation
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

# Include dependency files
-include $(OBJ:.o=.d)

# Debug build
debug:
	$(MAKE) CXXFLAGS="-std=c++17 -g -O0 -Wall -Wextra"

# Clean
clean:
	rm -f $(TARGET) $(OBJ) $(OBJ:.o=.d)

.PHONY: all clean debug

