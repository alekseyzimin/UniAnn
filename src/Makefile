# Compiler and flags
CXX      := g++
CXXFLAGS := -std=c++17 -O3 -Wall -Wextra -march=native
LDFLAGS  := 

# Target
U_TARGET := uniann 
U_SRC    := uniann.cpp
U_OBJ    := $(U_SRC:.cpp=.o)
C_TARGET := compute_markov_scores
C_SRC    := compute_markov_scores.cpp
C_OBJ    := $(C_SRC:.cpp=.o)


# Default rule
all: $(U_TARGET) $(C_TARGET)

# Link
$(U_TARGET): $(U_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
$(C_TARGET): $(C_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# Compile with automatic dependency generation
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

# Include dependency files
-include $(C_OBJ:.o=.d)  $(U_OBJ:.o=.d)

# Debug build
debug:
	$(MAKE) CXXFLAGS="-std=c++17 -g -O0 -Wall -Wextra"

# Clean
clean:
	rm -f $(U_TARGET) $(U_OBJ) $(U_OBJ:.o=.d) $(C_TARGET) $(C_OBJ) $(C_OBJ:.o=.d)

.PHONY: all clean debug

