#include "gmrt_newcorr.h"
#include "protocol.h"
#include "acqpsr.h"
#include "externalLibraries.h"
#include "frb_shm.h"
#ifndef CLASS_CORRELATOR
#define CLASS_CORRELATOR
class Correlator
{
public:
	void initializeReadSHM();
	void initializeReadSHM_SPOTLIGHT(char fileSHM);
	int initializeWriteSHM();
	int initializeWriteSHM_SPOTLIGHT(char fileSHM);
	void writeToSHM_SPOTLIGHT(unsigned char *rawData);
	void writeToSHM_SPOTLIGHT(unsigned char *rawData, struct timeval timestamp_gps);
	void writeToSHM(unsigned short int *rawData);
	void readFromSHM_SPOTLIGHT(unsigned char *rawData);
	void readFromSHM(unsigned short int *rawData);
	void writeToSHM(short int *rawData, char *header);
	void copyHeaderInfo();
	Correlator(int _nchan, float _sampling);
	Correlator(DasHdrType *_dataHdrRead, DataBufType *_dataBufferRead);

private:
	static DasHdrType *dataHdrWrite;
	static DataBufType *dataBufferWrite;
	static DataHeader *dataHdrWriteFRB;
	static DataBuffer *dataBufferWriteFRB;
	static DataTabType *dataTabWrite;

	static int recNumWrite;

	static DasHdrType *dataHdrRead;
	static DataBufType *dataBufferRead;

	static DataHeader *dataHdrReadFRB;
	static DataBuffer *dataBufferReadFRB;

	static DataTabType *dataTabRead;
	static int recNumRead;
	static long int currentReadBlock;
	int DataOff;
	char debug;
	int nchan;
	float sampling;
};
#endif