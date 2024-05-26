#include <AcquireData.h>

// Declaring static variables:
Information AcquireData::info;
long int AcquireData::eof;
DasHdrType *AcquireData::dataHdr;
DataBufType *AcquireData::dataBuffer;
DataTabType *AcquireData::dataTab;
GlobalInfoType *AcquireData::dataBufferProcess;
unsigned short *AcquireData::zeros;
int AcquireData::recNum = 0;
int AcquireData::remainingData = 0;
int AcquireData::currentReadBlock = 0;
long int AcquireData::curPos = 0;
double AcquireData::totalError = 0.0;
float AcquireData::buffSizeSec;
int AcquireData::nbuff;
long int AcquireData::bufferBlockLengthProcess;
struct timeval *AcquireData::startTimeStamp;
Correlator *AcquireData::shmInterface;

/*******************************************************************
 *CONSTRUCTOR: AcquireData::AcquireData(Information _info)
 *Information _info : contains all parameters.
 *Invoked on first creation of an object of this type.
 ********************************************************************/
AcquireData::AcquireData(Information _info)
{
	fprintf(stderr, "Inside AcquireData(Information _info) Constructor.\n");
	info = _info;
	blockIndex = 0;
	hasReachedEof = 0;
	totalError = 0.0;
	remainingData = 0;
	// finding length of data file or initializing SHM according to the mode of operation
	if (info.doReadFromFile)
	{
		fprintf(stderr, "Inside AcquireData(Information _info) Constructor, info.filepath is: %s\n", info.filepath);
		ifstream datafile;
		datafile.open(info.filepath + to_string(0), ios::binary);
		if (!datafile.is_open())
		{
			fprintf(stderr, "Raw data file not found!\n");
			exit(1);
		}
		datafile.seekg(0, ios::end);
		eof = datafile.tellg();
		datafile.close();
		if (eof == 0)
		{
			fprintf(stderr, "DATA FILE EMPTY!\n");
			exit(1);
		}
	}
	// else
	// initializeSHM();
}

/*******************************************************************
 *CONSTRUCTOR: AcquireData::AcquireData(int _blockIndex)
 *int _blockIndex : current window number.
 ********************************************************************/
AcquireData::AcquireData(int _blockIndex)
{
	blockIndex = _blockIndex;
	hasReachedEof = 0;
	rawData = NULL;
	rawDataFloat = NULL;
	rawDataPolar = NULL;
	splittedRawData = NULL;
	nBuffTaken = 0;
	if (info.isInline)
		headerInfo = new char[4096 * nbuff];
	headerInfo = NULL;
}

/*******************************************************************
 *DESTRUCTOR: AcquireData::~AcquireData()
 *Free all memory
 ********************************************************************/
AcquireData::~AcquireData()
{
	switch (info.sampleSizeBytes)
	{
	case 1:
		if (!info.doPolarMode)
		{
			delete[] rawDataChar;
		}
		else
			delete[] rawDataCharPolar;
		break;
	case 2:
		if (!info.doPolarMode)
			delete[] rawData;
		else
			delete[] rawDataPolar;
		break;
	case 4:
		delete[] rawDataFloat;
		break;
	}

	// delete[] splittedRawData;
}

/*******************************************************************
 *FUNCTION: void BasicAnalysis::readData()
 *Wrapper to decide which method to call to acquire data
 *******************************************************************/
void AcquireData::readData(int iBeam)
{
	fprintf(stderr, "AcquireData::readData(iBeam: %d)\n", iBeam);
	double er = info.blockSizeSec / info.samplingInterval;
	er = er - (long)er;
	totalError += er;
	blockLength = info.blockSizeSamples * info.timeIntFactor;
	nChannel = info.noOfChannels * info.freqIntFactor;
	if (!info.isInline && info.timeIntFactor == 1 && (info.shmID <= 3) && totalError >= 1.0)
	{
		blockLength++;
		totalError--;
	}
	if (info.doReadFromFile)
		readDataFromMultipleFiles(iBeam);
	else
	{
		if (info.shmID == 4 or info.shmID == 5)
		{
			long samplesToTake = nChannel * info.noOfPol * blockLength * info.sampleSizeBytes;
			rawDataChar = new unsigned char[blockLength * info.noOfPol * nChannel];
			shmInterface->readFromSHM_SPOTLIGHT(rawDataChar);
		}
		else if (info.shmID == 7)
			readFromSHM_SPOTLIGHT();
	}
}

/*******************************************************************
 *FUNCTION: void AcquireData::initializeSHM()
 *Initializes shared memory of GSB/GWB
 *******************************************************************/
void AcquireData::initializeSHM()
{
	fprintf(stderr, "AcquireData::initializeSHM()\n");
	int iddataHdr, idDataBuffer;
	if (info.isInline && !info.doFRB && info.shmID <= 3)
	{
		buffSizeSec = (info.samplingInterval * info.timeIntFactor * info.freqIntFactor * (dataBuffer->blocksize - 4096) / (info.sampleSizeBytes * info.noOfChannels));
		nbuff = (int)floor(info.blockSizeSec / buffSizeSec);
		info.blockSizeSec = nbuff * buffSizeSec;
		info.blockSizeSamples = (int)round(info.blockSizeSec / info.samplingInterval);
	}
	else if (info.doFRB || info.shmID >= 4)
	{
		buffSizeSec = info.samplingInterval * info.timeIntFactor * info.freqIntFactor * ((long)DataSize / (info.noOfChannels));
		nbuff = 1;
		info.blockSizeSec = nbuff * buffSizeSec;
		info.blockSizeSamples = info.timeIntFactor * info.freqIntFactor * (long)DataSize / info.noOfChannels;
	}
	switch (info.shmID)
	{
	case 4: // FRB shm
		// file shm
		fprintf(stderr, "Attempting to attach to file simulator FRB SHM.\n");
		shmInterface = new Correlator(info.noOfChannels, info.samplingInterval);
		shmInterface->initializeReadSHM_SPOTLIGHT(1);
		break;
	case 5: // FRB file shm
		// file shm
		fprintf(stderr, "Attempting to attach to FRB SHM.\n");
		shmInterface = new Correlator(info.noOfChannels, info.samplingInterval);
		shmInterface->initializeReadSHM_SPOTLIGHT(0);
		break;
	case 7:
		// process psr shm
		fprintf(stderr, "Attempting to attach to SPOTLIGHT SHM.\n");
		initializeSHM_SPOTLIGHT();
		break;
	}
}

void AcquireData::initializeSHM_SPOTLIGHT()
{
	fprintf(stderr, "AcquireData::initializeSHM_SPOTLIGHT()\n");
	bufferBlockLengthProcess = TEL_SHM_BlockSize; // 3.175 MB for 800 time samples and 4096 channel; 1.31072 ms.

	int local_ExtraBuf = 64;
	int CurrShmDataSize = (((MaxRecs + 1) * bufferBlockLengthProcess + local_ExtraBuf) / PageSize + 1) * PageSize;
	int local_ShmDataOff = ((sizeof(GlobalInfoType) + local_ExtraBuf) / PageSize + 1) * PageSize;
	int CurrShmSize = local_ShmDataOff + CurrShmDataSize;
	int idDataBuffer = shmget(ShmKey, CurrShmSize, 0);
	fprintf(stderr, "idDataBuffer: %d\n", idDataBuffer);
	if (idDataBuffer < 0)
	{
		exit(1);
	}

	dataBufferProcess = (GlobalInfoType *)shmat(idDataBuffer, 0, SHM_RDONLY);
	if ((void *)dataBufferProcess == (void *)-1)
	{
		fprintf(stderr, "Failed to attach to SPOTLIGHT SHM.\n");
		exit(1);
	}

	cout << "Attached to SPOTLIGHT shared memory: " << dataBufferProcess << endl;
	fprintf(stderr, "Max no of blocks: %d", MaxRecs);

	recNum = (dataBufferProcess->rec_ind + MaxRecs - 3) % MaxRecs;
	currentReadBlock = dataBufferProcess->rec[recNum].rec_seq;
}

int AcquireData::readFromSHM_SPOTLIGHT()
{
	fprintf(stderr, "AcquireData::readFromSHM_SPOTLIGHT()\n");
	int DataOff = 4096;

	long samplesToTake = nChannel * info.noOfPol * blockLength * info.sampleSizeBytes;
	rawDataChar = new unsigned char[blockLength * info.noOfPol * nChannel];
	long int fetched = 0;
	ofstream meanFile;
	meanFile.open("realTimeWarning.gpt", ios::app);
	ofstream recordFile;
	recordFile.open("blockLossRecord.gpt", ios::app);
	char *shmp = (char *)dataBufferProcess;
	while (fetched < samplesToTake)
	{
		timeWaitTime += omp_get_wtime();
		int flag = 0;

		while (dataBufferProcess->rec[recNum].rec_flag & Marked)
		{
			usleep(2000);
			if (flag == 0)
			{
				fprintf(stderr, "Waiting for flag...");
				flag = 1;
			}
		}
		if (flag == 1)
			fprintf(stderr, "Got a flag. Ready!");

		timeWaitTime -= omp_get_wtime();
		currentReadBlock = dataBufferProcess->rec[recNum].rec_seq;

		int ind = dataBufferProcess->blk_ind > 0 ? dataBufferProcess->blk_ind - 1 : MaxBLK - 1;

		int curInd = dataBufferProcess->rec_ind - 1;
		if (curInd < 0)
			curInd = 0;
		if (dataBufferProcess->rec[curInd].rec_seq >= currentReadBlock + MaxRecs - 2)
		{
			meanFile << "Lag in block " << blockIndex << endl;
			meanFile << "recNum = " << recNum << ", Reading Sequence: " << currentReadBlock << "dataBufferProcess->blk_seq[ind]" << dataBufferProcess->blk_seq[ind] << endl;
			fprintf(stderr, "recNum: %d; currentReadBlock: %d; rec_seq: %u\n", recNum, currentReadBlock, dataBufferProcess->rec[curInd].rec_seq);
			fprintf(stderr, "Difference, i.e., lag: %u\n", dataBufferProcess->rec[curInd].rec_seq - currentReadBlock);
			fprintf(stderr, "Max no of blocks - 2 = %d\n", MaxRecs - 2);
			fprintf(stderr, "Processing lagged behind...\n");

			recNum = (dataBufferProcess->rec_ind + MaxRecs - 2) % MaxRecs;
			currentReadBlock = dataBufferProcess->rec[recNum].rec_seq;
		}
		timestamp_gps = dataBufferProcess->rec[recNum].timestamp_gps;
		fprintf(stderr, "recNum: %d; currentReadBlock: %d\n", recNum, currentReadBlock); // Collect's Sequence: dataBuffer->cur_block-1;
		unsigned char *bufptr = (unsigned char *)(shmp + dataBufferProcess->rec[recNum].beg_off + remainingData);
		if (samplesToTake - fetched >= bufferBlockLengthProcess - remainingData)
		{
			if (info.isInline)
			{
				// memcpy(headerInfo+nBuffTaken*DataOff, dataBuffer->buf+dataTab[recNum].rec*(dataBuffer->blocksize), DataOff);
				nBuffTaken++;
			}
			memcpy(rawDataChar + fetched / sizeof(char), bufptr, bufferBlockLengthProcess - remainingData);

			fetched += (bufferBlockLengthProcess - remainingData);
			recNum = (recNum + 1) % MaxRecs;
			remainingData = 0;
		}
		else
		{
			memcpy(rawDataChar + fetched / sizeof(char), bufptr, samplesToTake - fetched);
			remainingData = (samplesToTake - fetched);
			break;
		}
	}
	curPos += samplesToTake;
	meanFile.close();
	return 1;
}

/*******************************************************************
 *FUNCTION: void AcquireData::readDataFromFile()
 *In offline mode reads raw data from a file
 *******************************************************************/
void AcquireData::readDataFromMultipleFiles(int iBeam)
{
	fprintf(stderr, "Inside AcquireData::readDataFromMultipleFiles(iBeam: %d\n)", iBeam);
	int c, i;
	long int blockSizeBytes = nChannel * info.noOfPol * blockLength * info.sampleSizeBytes; // Number of bytes to read in each block	// Number of bytes that have already been read

	ifstream datafile;
	fprintf(stderr, "File to be read, info.filepath + to_string(iBeam), is: ");
	cout << info.filepath + to_string(iBeam) << endl;
	datafile.open(info.filepath + to_string(iBeam), ios::binary);
	datafile.seekg(curPos, ios::beg);
	// logic to handle reading last block
	if (curPos + blockSizeBytes > eof)
	{
		blockSizeBytes = eof - curPos;
		blockLength = blockSizeBytes / (info.sampleSizeBytes * nChannel * info.noOfPol); // The number of time samples
		hasReachedEof = 1;
	}

	rawDataChar = new unsigned char[blockSizeBytes / info.sampleSizeBytes];
	datafile.read((char *)rawDataChar, blockSizeBytes);

	datafile.close();
	curPos += blockSizeBytes;
}

float u16tofloat(short x)
{
	fprintf(stderr, "Inside AcquireData::u16tofloat(short x)\n");
	union
	{
		float f;
		int i;
	} u;
	u.f = 0.00f;
	u.i |= x;
	return u.f;
}

/**********************************************************************
*FUNCTION: void AcquireData::splitRawData()
*If on polar mode, splitting of raw data into four polarization channels
*is done here. GMRT polarization data is read out in the following fomat:
T1_C1_P T1_C1_Q T1_C1_R T1_C1_S T1_C2_P T1_C2_Q T1_C2_R T1_C2_S...
T1_CN_S T2_C1_P ...........TN_CN_S
where P,Q,R,S are the four polarizations, Cn denotes the nth channel
and Tn the nth time sample. So T15_C25_Q means the Q-polarization data
of channel 25 of the 15th time sample.
**********************************************************************/
void AcquireData::splitRawData()
{
	fprintf(stderr, "Inside AcquireData::splitRawData()\n");
	splittedRawData = new float *[info.noOfPol];
	float **ptrSplittedRawData = new float *[info.noOfPol];
	int nChan = info.noOfChannels * info.freqIntFactor;
	long int length = blockLength * nChan;
	for (int k = 0; k < info.noOfPol; k++)
	{
		splittedRawData[k] = new float[length];
		ptrSplittedRawData[k] = splittedRawData[k];
	}
	/*******************************************************************
	 *Handles different data types. GMRT data is mostly of type short while
	 *certain processed data maybe floating point.
	 *******************************************************************/

	float *ptrRawDataFloat = rawDataFloat;
	switch (info.sampleSizeBytes)
	{
	case 1:
		if (!info.doPolarMode)
		{
			unsigned char *ptrRawData = rawDataChar;
			for (int j = 0; j < blockLength; j++) // Refer to the GMRT polarization data format in Function description.
			{
				for (int i = 0; i < nChan; i++, ptrRawData++, ptrSplittedRawData[0]++)
				{
					*(ptrSplittedRawData[0]) = (*ptrRawData);
				}
			}
		}
		else
		{
			char *ptrRawData = rawDataCharPolar;
			for (int j = 0; j < blockLength; j++) // Refer to the GMRT polarization data format in Function description.
			{
				for (int i = 0; i < nChan; i++)
				{
					*(ptrSplittedRawData[0]++) = (*ptrRawData++);
					*(ptrSplittedRawData[1]++) = (*ptrRawData++);
					*(ptrSplittedRawData[2]++) = (*ptrRawData++);
					*(ptrSplittedRawData[3]++) = (*ptrRawData++);
				}
			}
		}
		break;
	case 2:
		if (!info.doPolarMode)
		{
			unsigned short int *ptrRawData = rawData;
			for (int j = 0; j < blockLength; j++) // Refer to the GMRT polarization data format in Function description.
			{
				for (int i = 0; i < nChan; i++, ptrRawData++, ptrSplittedRawData[0]++)
				{
					*(ptrSplittedRawData[0]) = (*ptrRawData);
				}
			}
		}
		else
		{
			short int *ptrRawData = rawDataPolar;
			for (int j = 0; j < blockLength; j++) // Refer to the GMRT polarization data format in Function description.
			{
				for (int i = 0; i < nChan; i++)
				{
					*(ptrSplittedRawData[0]++) = (*ptrRawData++);
					*(ptrSplittedRawData[1]++) = (*ptrRawData++);
					*(ptrSplittedRawData[2]++) = (*ptrRawData++);
					*(ptrSplittedRawData[3]++) = (*ptrRawData++);
				}
			}
		}
		break;
	case 4:
		for (int j = 0; j < blockLength; j++) // Refer to the GMRT polarization data format in Function description.
		{
			for (int i = 0; i < nChan; i++)
			{
				for (int k = 0; k < info.noOfPol; k++, ptrRawDataFloat++)
				{
					(*ptrSplittedRawData[k]) = (*ptrRawDataFloat);
					ptrSplittedRawData[k]++;
				}
			}
		}
		break;
	}
	delete[] ptrSplittedRawData;
}

void AcquireData::averageRawData(int nIntTime, int nIntFreq)
{
	fprintf(stderr, "Inside AcquireData::averageRawData(int nIntTime, int nIntFreq)\n");
	long int blockSizeBytes = nChannel * info.noOfPol * blockLength * info.sampleSizeBytes;
	float **splittedRawDataAvg = new float *[info.noOfPol];
	float *ptrSplittedRawDataAvg;
	float *ptrSplittedRawData;
	long int length = blockLength / info.timeIntFactor * info.noOfChannels;
	for (int k = 0; k < info.noOfPol; k++)
	{
		splittedRawDataAvg[k] = new float[length];
	}
	int newNChan = info.noOfChannels;
	int newNlen = blockLength / nIntTime;
	for (int k = 0; k < info.noOfPol; k++)
	{
		ptrSplittedRawDataAvg = splittedRawDataAvg[k];
		for (int i = 0; i < newNlen * newNChan; i++, ptrSplittedRawDataAvg++)
		{
			*ptrSplittedRawDataAvg = 0;
		}
		ptrSplittedRawDataAvg = splittedRawDataAvg[k];
		ptrSplittedRawData = splittedRawData[k];
		for (int j = 0; j < newNlen; j++)
		{
			for (int n = 0; n < nIntTime; n++)
			{
				ptrSplittedRawDataAvg = &splittedRawDataAvg[k][j * newNChan];
				for (int i = 0; i < newNChan; i++, ptrSplittedRawDataAvg++)
				{
					for (int m = 0; m < nIntFreq; m++, ptrSplittedRawData++)
					{
						*ptrSplittedRawDataAvg += (*ptrSplittedRawData) / (nIntTime * nIntFreq);
					}
				}
			}
		}
		delete[] splittedRawData[k];
	}
	splittedRawData = splittedRawDataAvg;
}