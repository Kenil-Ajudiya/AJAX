#ifndef SOFTCORR_H
#define SOFTCORR_H

#define DasHeaderKey 2031
#define DasBufferKey 2032

#define DasHeaderKey_SIM 3031
#define DasBufferKey_SIM 3032

// #include <sys/ipc.h>
// #include <sys/shm.h>
// #include <sys/types.h>
#define process_psr_SHM_BlockSize (512 * 4096)
// 32 blocks of the multi-beam version of TEL_SHM (process_psr), each having 800 time samples and 4096 channels; 8-bit integers.
// Time resolution at this stage is 1.31072 ms, which makes one block of the multi-beam version of TEL_SHM (process_psr) 1.048576 s long.
#define TEL_SHM_BlockSize (800 * 4096)
#define DataSize (32 * TEL_SHM_BlockSize) // Size of one block of FRB_BEAM_SHM = 100 MB; 33.554432 s.
#define MaxDataBlocks 12 // There will be 12 blocks, each of the above mentioned size, totalling 1.171875 GB; 402.653184 s.
#define TimeSize sizeof(double)

typedef struct
{
  unsigned int flag, curBlock, curRecord, blockSize;
  int overFlow;
  double pcTime[MaxDataBlocks], dataTime[MaxDataBlocks];
  unsigned char data[(long)(DataSize) * (long)(MaxDataBlocks)];
} DataBuffer;

typedef struct
{
  unsigned int active, status;
  double pcTime, dataTime, refTime;
  struct timeval timestamp[MaxDataBlocks];
  struct timeval timestamp_gps[MaxDataBlocks];
  double blk_nano[MaxDataBlocks];
} DataHeader;

// DataHeader *dataHdr;
// DataBuffer *dataBuf;

#endif // SOFTCORR_H