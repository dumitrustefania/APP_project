CXX = g++
CXXFLAGS = -std=c++17
TARGET_SEQUENTIAL = sequential
TARGET_OPENMP = openmp
TARGET_MPI = mpi
TARGET_HYBRID = hybrid

INPUT_DIR = inputs
OUTPUT_DIR = outputs
OUTPUT_DIR_SEQUENTIAL = $(OUTPUT_DIR)/sequential
OUTPUT_DIR_OPENMP = $(OUTPUT_DIR)/openmp
OUTPUT_DIR_MPI = $(OUTPUT_DIR)/mpi
OUTPUT_DIR_HYBRID = $(OUTPUT_DIR)/hybrid

all: $(TARGET_SEQUENTIAL) $(TARGET_OPENMP) $(TARGET_MPI) $(TARGET_HYBRID)

$(TARGET_SEQUENTIAL): src/sequential.cpp
	$(CXX) $(CXXFLAGS) -o $(TARGET_SEQUENTIAL) src/sequential.cpp

$(TARGET_OPENMP): src/openmp.cpp
	$(CXX) $(CXXFLAGS) -fopenmp -o $(TARGET_OPENMP) src/openmp.cpp

$(TARGET_MPI): src/mpi.cpp
	mpicxx $(CXXFLAGS) -o $(TARGET_MPI) src/mpi.cpp

$(TARGET_HYBRID): src/hybrid.cpp
	mpicxx $(CXXFLAGS) -fopenmp -o $(TARGET_HYBRID) src/hybrid.cpp

run_sequential: $(TARGET_SEQUENTIAL)
	@mkdir -p $(OUTPUT_DIR_SEQUENTIAL)
	./$(TARGET_SEQUENTIAL) $(INPUT_DIR)/small.in
	./$(TARGET_SEQUENTIAL) $(INPUT_DIR)/medium.in
	./$(TARGET_SEQUENTIAL) $(INPUT_DIR)/large.in

run_openmp: $(TARGET_OPENMP)
	@mkdir -p $(OUTPUT_DIR_OPENMP)
	./$(TARGET_OPENMP) $(INPUT_DIR)/small.in 4
	./$(TARGET_OPENMP) $(INPUT_DIR)/medium.in 4
	./$(TARGET_OPENMP) $(INPUT_DIR)/large.in 4

run_mpi: $(TARGET_MPI)
	@mkdir -p $(OUTPUT_DIR_MPI)
	mpirun -np 4 ./$(TARGET_MPI) $(INPUT_DIR)/small.in
	mpirun -np 4 ./$(TARGET_MPI) $(INPUT_DIR)/medium.in
	mpirun -np 4 ./$(TARGET_MPI) $(INPUT_DIR)/large.in

run_hybrid: $(TARGET_HYBRID)
	@mkdir -p $(OUTPUT_DIR_HYBRID)
	mpirun -np 2 ./$(TARGET_HYBRID) $(INPUT_DIR)/small.in 4
	mpirun -np 2 ./$(TARGET_HYBRID) $(INPUT_DIR)/medium.in 4
	mpirun -np 2 ./$(TARGET_HYBRID) $(INPUT_DIR)/large.in 4

run: run_sequential run_openmp run_mpi run_hybrid

clean:
	@rm -f $(TARGET_SEQUENTIAL)
	@rm -rf $(OUTPUT_DIR_SEQUENTIAL)
	@rm -f $(TARGET_OPENMP)
	@rm -rf $(OUTPUT_DIR_OPENMP)
	@rm -f $(TARGET_MPI)
	@rm -rf $(OUTPUT_DIR_MPI)
	@rm -f $(TARGET_HYBRID)
	@rm -rf $(OUTPUT_DIR_HYBRID)
	@rm -rf $(OUTPUT_DIR)

.PHONY: all run_sequential run_openmp run_mpi run_hybrid clean
