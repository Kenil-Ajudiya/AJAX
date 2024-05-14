using namespace std;
#include "externalLibraries.h"
#include "Information.h"
#include "SHM.h"

#ifndef CLASS_ACQUIREDATA
#define CLASS_ACQUIREDATA
/*******************************************************************
 *CLASS:	AcquireData
 *Contains functions to read data either from a file or SHM (in case
 *of real time operation). It also has a functions to split read data
 *into four polarizations.
 *******************************************************************/
class AcquireData
{
public:
	// Static variables:
	static Information info;
	static long int eof;	  // Length of file in bytes
	static double totalError; // Cumulative uncorrected difference between true window width and required width
	static long int curPos;	  // Total number of samples read.
	// Static SHM Variables:
	static DasHdrType *dataHdr;
	static DataBufType *dataBuffer;
	static DataTabType *dataTab;
	static GlobalInfoType *dataBufferProcess;
	static unsigned short *zeros;
	static int recNum;
	static int currentReadBlock;
	static int remainingData;
	static timeval *startTimeStamp;
	static float buffSizeSec;
	static int nbuff;
	static Correlator *shmInterface;
	static long int bufferBlockLengthProcess;
	// Variables:
	unsigned char *rawDataChar; // 1-byte integer read data is stored here
	char *rawDataCharPolar;		// 1-byte integer read data is stored here
	unsigned short *rawData;	// 2-byte integer read data is stored here
	short *rawDataPolar;		// 2-byte integer read data is stored here
	float *rawDataFloat;		// 4-byte floating point read data
	float **splittedRawData;	// Pointers to 4 polarization data in case of polar modeIndex. In normal run it just has one pointer to rawData.
	char hasReachedEof;			// Flag to let the program know when end of file has been reached
	long int blockLength;		// Length(in samples) of the current block
	int nChannel;				// Number of input channels
	int blockIndex;				// Current window number. (starting from 0)
	char *headerInfo;
	int nBuffTaken;
	double timeWaitTime; // Benchmark
	struct timeval timestamp_gps;

	// Functions:
	AcquireData(Information _info);			   // Constructor for first time intialization, the static variables are initialized here.
	AcquireData(int _blockIndex);			   // Constructor to intialize non-static me
	~AcquireData();							   // *****WHAT IS THIS FOR???*****
	void readData(int iBeam);				   // Calls appropiate data reading function for the multi-beam version.
	void readDataFromMultipleFiles(int iBeam); // Reads from  multiple files for multi-beam version.
	void initializeSHM();					   // Attaches to SHM
	void initializeSHM_SPOTLIGHT();			   // attaches to SPOTLIGHT SHM
	int readFromSHM_SPOTLIGHT();			   // Reads from process SHM
	void splitRawData();
	void averageRawData(int nIntTime, int nIntFreq);
};
#endif