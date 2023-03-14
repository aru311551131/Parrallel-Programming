CXX := mpicxx
RUN := mpirun
CXXFLAG := -O3 -g
SRC := $(notdir $(wildcard *.cc))
OBJ := $(patsubst %.cc,%.o,$(SRC))
TARGET := matmul

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(OBJ) -o $(TARGET) $(CXXFLAG)

%.o: %.cc
	$(CXX) -c $< -o $@ $(CXXFLAG)

clean:
	rm -f $(OBJ) $(TARGET)