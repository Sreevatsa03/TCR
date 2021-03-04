#include <wb.h>

#define wbCheck(stmt)                                                     \
  do {                                                                    \
    cudaError_t err = stmt;                                               \
    if (err != cudaSuccess) {                                             \
      wbLog(ERROR, "Failed to run stmt ", #stmt);                         \
      wbLog(ERROR, "Got CUDA error ...  ", cudaGetErrorString(err));      \
      return -1;                                                          \
    }                                                                     \
  } while (0)

#define TILE_WIDTH 16

// Compute C = A * B
__global__ void matrixMultiplyShared(float *A, float *B, float *C,
                                     int numARows, int numAColumns,
                                     int numBRows, int numBColumns,
                                     int numCRows, int numCColumns) {
  //@@ Insert code to implement matrix multiplication here
  //@@ You have to use tiling with shared memory for arbitrary size
  float CValue = 0;

  int row = blockIdx.y * TILE_WIDTH + threadIdx.y;
  int col = blockIdx.x * TILE_WIDTH + threadIdx.x;

  __shared__ float As[TILE_WIDTH][TILE_WIDTH];
  __shared__ float Bs[TILE_WIDTH][TILE_WIDTH];

  for (int i = 0; i < (TILE_WIDTH + numAColumns - 1)/TILE_WIDTH; i++) {

    if (i * TILE_WIDTH + threadIdx.x < numAColumns && row < numARows)
      As[threadIdx.y][threadIdx.x] = A[row * numAColumns + i*TILE_WIDTH + threadIdx.x];
    else
      As[threadIdx.y][threadIdx.x] = 0.0;

    if (i * TILE_WIDTH + threadIdx.y < numBRows && col < numBColumns)
      Bs[threadIdx.y][threadIdx.x] = B[(i * TILE_WIDTH + threadIdx.y) * numBColumns + col];
    else
      Bs[threadIdx.y][threadIdx.x] = 0.0;

    __syncthreads();

    for (int j = 0; j < TILE_WIDTH; ++j)
      CValue += As[threadIdx.y][j] * Bs[j][threadIdx.x];

    __syncthreads();
  }

  if (row < numCRows && col < numCColumns)
    C[((blockIdx.y * blockDim.y + threadIdx.y) * numCColumns) + (blockIdx.x * blockDim.x) + threadIdx.x] = CValue;
}

int main(int argc, char **argv) {
  wbArg_t args;
  float *hostA; // The A matrix
  float *hostB; // The B matrix
  float *hostC; // The output C matrix
  float *deviceA; // A matrix on device
  float *deviceB; // B matrix on device
  float *deviceC; // C matrix on device
  int numARows;    // number of rows in the matrix A
  int numAColumns; // number of columns in the matrix A
  int numBRows;    // number of rows in the matrix B
  int numBColumns; // number of columns in the matrix B
  int numCRows;    // number of rows in the matrix C(you have to set this)
  int numCColumns; // number of columns in the matrix C (you have to set
                   // this)

  args = wbArg_read(argc, argv);

  wbTime_start(Generic, "Importing data and creating memory on host");
  hostA = (float *)wbImport(wbArg_getInputFile(args, 0), &numARows,
                            &numAColumns);
  hostB = (float *)wbImport(wbArg_getInputFile(args, 1), &numBRows,
                            &numBColumns);
                            
  //@@ Set numCRows and numCColumns
  numCRows    = numARows;   // set to correct value
  numCColumns = numBColumns;   // set to correct value
  //@@ Allocate the hostC matrix
  hostC = (float *) malloc(sizeof(float) * numCRows * numCColumns);
  
  wbTime_stop(Generic, "Importing data and creating memory on host");

  wbLog(TRACE, "The dimensions of A are ", numARows, " x ", numAColumns);
  wbLog(TRACE, "The dimensions of B are ", numBRows, " x ", numBColumns);
  wbLog(TRACE, "The dimensions of C are ", numCRows, " x ", numCColumns);
  
  wbTime_start(GPU, "Allocating GPU memory.");
  //@@ Allocate GPU memory here
  cudaMalloc((void **) &deviceA, sizeof(float)*numARows*numAColumns);
  cudaMalloc((void **) &deviceB, sizeof(float)*numBRows*numBColumns);
  cudaMalloc((void **) &deviceC, sizeof(float)*numCRows*numCColumns);

  wbTime_stop(GPU, "Allocating GPU memory.");

  wbTime_start(GPU, "Copying input memory to the GPU.");
  //@@ Copy memory to the GPU here
  cudaMemcpy(deviceA, hostA, sizeof(float)*numARows*numAColumns, cudaMemcpyHostToDevice);
  cudaMemcpy(deviceB, hostB, sizeof(float)*numBRows*numBColumns, cudaMemcpyHostToDevice);

  wbTime_stop(GPU, "Copying input memory to the GPU.");

  //@@ Initialize the grid and block dimensions here
  // note that TILE_WIDTH is set to 16 on line number 13. 
  dim3 dimBlock(TILE_WIDTH, TILE_WIDTH, 1);    
  dim3 dimGrid((numCColumns/TILE_WIDTH) + 1, (numCRows/TILE_WIDTH) + 1, 1);
  
  wbTime_start(Compute, "Performing CUDA computation");
  //@@ Launch the GPU Kernel here
  matrixMultiplyShared<<<dimGrid, dimBlock>>>(deviceA, deviceB, deviceC, numARows, numAColumns, numBRows, numBColumns, numCRows, numCColumns);

  cudaDeviceSynchronize();
  wbTime_stop(Compute, "Performing CUDA computation");

  wbTime_start(Copy, "Copying output memory to the CPU");
  //@@ Copy the GPU memory back to the CPU here
  cudaMemcpy(hostC, deviceC, sizeof(float)*numCRows*numCColumns, cudaMemcpyDeviceToHost);

  wbTime_stop(Copy, "Copying output memory to the CPU");

  wbTime_start(GPU, "Freeing GPU Memory");
  //@@ Free the GPU memory here
  cudaFree(deviceA);
  cudaFree(deviceB);
  cudaFree(deviceC);
  
  wbTime_stop(GPU, "Freeing GPU Memory");

  wbSolution(args, hostC, numCRows, numCColumns);

  free(hostA);
  free(hostB);
  free(hostC);

  return 0;
}
