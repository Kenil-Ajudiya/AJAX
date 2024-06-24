/*************************************************************
shm library
This library is aimed at interfacing between various
kinds of shared memory structures used by different
pulsar codes at the GMRT.
				-Aditya Chowdhury, 14th March 2016
************************************************************/
#include "SHM.h"
using namespace std;

DasHdrType *Correlator::dataHdrWrite;
DataBufType *Correlator::dataBufferWrite;
DataTabType *Correlator::dataTabWrite;
DataHeader *Correlator::dataHdrWriteFRB;
DataBuffer *Correlator::dataBufferWriteFRB;
DasHdrType *Correlator::dataHdrRead;
DataBufType *Correlator::dataBufferRead;
DataTabType *Correlator::dataTabRead;
DataHeader *Correlator::dataHdrReadFRB;
DataBuffer *Correlator::dataBufferReadFRB;
int Correlator::recNumRead = 0;
int Correlator::recNumWrite = 0;
long int Correlator::currentReadBlock = 0;

Correlator::Correlator(DasHdrType *_dataHdrRead, DataBufType *_dataBufferRead)
{
	DataOff = 4096;
	debug = 1;
	dataHdrRead = _dataHdrRead;
	dataBufferRead = _dataBufferRead;
}

Correlator::Correlator(int _nchan, float _sampling)
{
	DataOff = 4096;
	debug = 1;
	nchan = _nchan;
	sampling = _sampling;
	fprintf(stderr, "\nOutput nchan: %d; sampling: %d\n", nchan, sampling);
	dataHdrRead = NULL;
	dataBufferRead = NULL;
}

/*******************************************************************
 *FUNCTION: void AcquireData::initializeSHM()
 *Initializes collect_psr shared memory of GSB/GWB
 *******************************************************************/
void Correlator::initializeReadSHM()
{
	fprintf(stderr, "Inside Correlator::initializeReadSHM()\n");
	int idDataHdr, idDataBuffer;
	idDataHdr = shmget(DAS_H_KEY, sizeof(DasHdrType), 644);
	idDataBuffer = shmget(DAS_D_KEY, sizeof(DataBufType), 0);
	fprintf(stderr, "idDataHdr: %d\nidDataBuffer: %d\n", idDataHdr, idDataBuffer);

	if (idDataHdr < 0 || idDataHdr < 0)
	{
		exit(1);
	}
	dataHdrRead = (DasHdrType *)shmat(idDataHdr, 0, SHM_RDONLY);
	dataBufferRead = (DataBufType *)shmat(idDataBuffer, 0, SHM_RDONLY);
	if ((dataBufferRead) == (DataBufType *)-1)
	{
		fprintf(stderr, "Cannot attach to SHM!\n");
		exit(1);
	}
	cout << "Attached to shared memory:" << dataHdrRead << "," << dataBufferRead << endl;
	dataTabRead = dataBufferRead->dtab;
	fprintf(stderr, "Max no. of blocks: %d\n", dataBufferRead->maxblocks);

	// find a block, two blocks before the current block of the shm for reading data
	// if(dataBuffer->cur_rec > (dataBuffer->maxblocks)/2)
	recNumRead = (dataBufferRead->cur_rec - 2 + MaxDataBuf) % MaxDataBuf;
	currentReadBlock = dataTabRead[recNumRead].seqnum;
	fprintf(stderr, "Exiting Correlator::initializeReadSHM()\n");
}

/*******************************************************************
 *FUNCTION: void AcquireData::initializeSHM()
 *Initializes collect_psr shared memory of GSB/GWB
 *******************************************************************/
void Correlator::initializeReadSHM_SPOTLIGHT(char fileSHM)
{
	fprintf(stderr, "Inside Correlator::initializeReadSHM_SPOTLIGHT(char fileSHM)\n");
	int idDataHdr, idDataBuffer;
	if (!fileSHM)
	{
		idDataHdr = shmget(DasHeaderKey, sizeof(DataHeader), 644);
		idDataBuffer = shmget(DasBufferKey, sizeof(DataBuffer), 644);
	}
	else
	{
		idDataHdr = shmget(DasHeaderKey_SIM, sizeof(DataHeader), 644);
		idDataBuffer = shmget(DasBufferKey_SIM, sizeof(DataBuffer), 644);
	}
	fprintf(stderr, "idDataHdr: %d\nidDataBuffer: %d\n", idDataHdr, idDataBuffer);

	if (idDataHdr < 0 || idDataHdr < 0)
	{
		exit(1);
	}
	dataHdrReadFRB = (DataHeader *)shmat(idDataHdr, 0, SHM_RDONLY);
	dataBufferReadFRB = (DataBuffer *)shmat(idDataBuffer, 0, SHM_RDONLY);
	cout << "Attached to shared memory:" << dataHdrReadFRB << "," << dataBufferReadFRB << endl;
	fprintf(stderr, "Max no. of blocks: %d\n", MaxDataBlocks);

	if ((dataBufferReadFRB) == (DataBuffer *)-1)
	{
		fprintf(stderr, "Cannot attach to SHM!\n");
		exit(1);
	}
	/* find a block, two blocks before the current block of the shm for reading data */
	// if(dataBuffer->cur_rec > (dataBuffer->maxblocks)/2)
	recNumRead = (dataBufferReadFRB->curRecord - 2 + MaxDataBlocks) % MaxDataBlocks;
	currentReadBlock = dataBufferReadFRB->curBlock - 2;

	fprintf(stderr, "Exiting Correlator::initializeReadSHM_SPOTLIGHT(char fileSHM)\n");
}

void Correlator::copyHeaderInfo()
{
	fprintf(stderr, "Inside Correlator::copyHeaderInfo()\n");
	dataHdrWrite->active = dataHdrRead->active;
	dataHdrWrite->status = dataHdrRead->status;
	dataHdrWrite->scan = dataHdrRead->scan;
	dataHdrWrite->scan_off = dataHdrRead->scan_off;
	dataHdrWrite->corr = dataHdrRead->corr;
	dataHdrWrite->model = dataHdrRead->model;
	dataHdrWrite->BeamHeader = dataHdrRead->BeamHeader;
	dataBufferWrite->blocksize = dataBufferRead->blocksize;
	dataBufferWrite->maxblocks = dataBufferRead->maxblocks;
	fprintf(stderr, "Exiting Correlator::copyHeaderInfo()\n");
}

/*******************************************************************
 *FUNCTION: void AcquireData::initializeSHM()
 *Initializes collect_psr shared memory of GSB/GWB
 *******************************************************************/
int Correlator::initializeWriteSHM()
{
	fprintf(stderr, "Inside Correlator::initializeWriteSHM()\n");
	int idDataHdr, idDataBuffer;
	if (dataHdrRead == NULL)
	{
		idDataHdr = shmget(DAS_H_KEY_AJAX, sizeof(DasHdrType), IPC_CREAT | 0666);
		idDataBuffer = shmget(DAS_D_KEY_AJAX, sizeof(DataBufType), IPC_CREAT | 0666);
	}
	else
	{
		idDataHdr = shmget(DAS_H_KEY_AJAX_INLINE, sizeof(DasHdrType), IPC_CREAT | 0666);
		idDataBuffer = shmget(DAS_D_KEY_AJAX_INLINE, sizeof(DataBufType), IPC_CREAT | 0666);
	}
	fprintf(stderr, "idDataHdr: %d\nidDataBuffer: %d\n", idDataHdr, idDataBuffer);

	if (idDataHdr < 0 || idDataHdr < 0)
	{
		fprintf(stderr, "Error creating SHM.\n");
		exit(1);
	}
	dataHdrWrite = (DasHdrType *)shmat(idDataHdr, 0, 0);
	dataBufferWrite = (DataBufType *)shmat(idDataBuffer, 0, 0);
	if ((dataBufferWrite) == (DataBufType *)-1)
	{
		fprintf(stderr, "Cannot attach to SHM!\n");
		exit(1);
	}
	cout << "Attached to write shared memory:" << dataHdrWrite << "," << dataBufferWrite << endl;

	if (dataHdrRead == NULL)
	{
		dataBufferWrite->blocksize = 2 * nchan * int(BLOCKTIME / sampling) + DataOff;
		dataBufferWrite->maxblocks = int(DAS_BUFSIZE / (dataBufferWrite->blocksize));
	}
	else
	{
		copyHeaderInfo();
		fprintf(stderr, "Header info copied.\n");
	}
	dataTabWrite = dataBufferWrite->dtab;
	dataBufferWrite->cur_rec = 0;
	dataBufferWrite->cur_block = 0;
	recNumWrite = (dataBufferWrite->cur_rec) % MaxDataBuf;
	dataTabWrite[recNumWrite].seqnum = 0;
	dataTabWrite[recNumWrite].rec = 0;
	for (int i = 0; i < MaxDataBuf; i++)
		dataTabWrite[i].flag = 0;
	fprintf(stderr, "Max no. of blocks: %d\n", dataBufferWrite->maxblocks);
	fprintf(stderr, "dataBufferWrite->blocksize: %d\n", dataBufferWrite->blocksize);
	dataHdrWrite->status = DAS_START;
	return dataBufferWrite->blocksize - DataOff;
	fprintf(stderr, "Exiting Correlator::initializeWriteSHM()\n");
}

/*******************************************************************
 *FUNCTION: void AcquireData::initializeSHM()
 *Initializes collect_psr shared memory of GSB/GWB
 *******************************************************************/
int Correlator::initializeWriteSHM_SPOTLIGHT(char fileSHM)
{
	fprintf(stderr, "Inside Correlator::initializeWriteSHM_SPOTLIGHT()\n");
	int idDataHdr, idDataBuffer;
	idDataBuffer = shmget(DasBufferKey, 50 * sizeof(DataBuffer), IPC_CREAT | 0666);	// Un-hardcode the number 50 with the nBeams.
	idDataHdr = shmget(DasHeaderKey, 50 * sizeof(DataHeader), IPC_CREAT | 0666);	// Un-hardcode the number 50 with the nBeams.

	fprintf(stderr, "idDataHdr: %d\nidDataBuffer: %d\n", idDataHdr, idDataBuffer);

	if (idDataHdr < 0 || idDataHdr < 0)
	{
		fprintf(stderr, "Error creating SHM!\n");
		exit(1);
	}
	dataHdrWriteFRB = (DataHeader *)shmat(idDataHdr, 0, 0);
	dataBufferWriteFRB = (DataBuffer *)shmat(idDataBuffer, 0, 0);
	if (dataBufferWriteFRB == (DataBuffer *)-1)
	{
		fprintf(stderr, "Cannot atttach to SHM!\n");
		exit(1);
	}
	cout << "Attached to write shared memory:" << dataHdrWriteFRB << "," << dataBufferWriteFRB << endl;

	dataBufferWriteFRB->curRecord = 0;
	dataBufferWriteFRB->curBlock = 0;
	recNumWrite = (dataBufferWriteFRB->curRecord) % MaxDataBlocks;

	fprintf(stderr, "Max no. of blocks: %d\n", MaxDataBlocks);
	fprintf(stderr, "DataSize: %d\n", DataSize);
	dataHdrWriteFRB->active = 1;
	return DataSize;
	fprintf(stderr, "Exiting Correlator::initializeWriteSHM_SPOTLIGHT()\n");
}

void Correlator::writeToSHM_SPOTLIGHT(unsigned char *rawData)
{
	fprintf(stderr, "Inside Correlator::writeToSHM_SPOTLIGHT(unsigned char *rawData)\n");
	fprintf(stderr, "Starting memcpy to SHM.\n");
	fprintf(stderr, "Writing Record: %u\n", dataBufferWriteFRB->curRecord);
	fprintf(stderr, "Writing Block: %u\n", dataBufferWriteFRB->curBlock);
	fprintf(stderr, "DataSize: %d\n", DataSize);
	memcpy(dataBufferWriteFRB->data + (long)DataSize * (long)recNumWrite, rawData, DataSize);
	fprintf(stderr, "Done memcpy to SHM.\n");

	dataBufferWriteFRB->curRecord = (recNumWrite + 1) % MaxDataBlocks;
	dataBufferWriteFRB->curBlock += 1;
	recNumWrite = (recNumWrite + 1) % MaxDataBlocks;
	fprintf(stderr, "Exiting Correlator::writeToSHM_SPOTLIGHT(unsigned char *rawData)\n");
}

void Correlator::writeToSHM_SPOTLIGHT(unsigned char *rawData, struct timeval timestamp_gps)
{
	fprintf(stderr, "Inside Correlator::writeToSHM_SPOTLIGHT(unsigned char *rawData, struct timeval timestamp_gps)\n");
	fprintf(stderr, "Starting memcpy to SHM.\n");
	fprintf(stderr, "Writing Record: %u\n", dataBufferWriteFRB->curRecord);
	fprintf(stderr, "Writing Block: %u\n", dataBufferWriteFRB->curBlock);
	fprintf(stderr, "DataSize: %d\n", DataSize);
	memcpy(dataBufferWriteFRB->data + (long)DataSize * (long)recNumWrite, rawData, DataSize);
	fprintf(stderr, "Done memcpy to SHM.\n");

	dataHdrWriteFRB->timestamp_gps[recNumWrite] = timestamp_gps;
	dataBufferWriteFRB->curRecord = (recNumWrite + 1) % MaxDataBlocks;
	dataBufferWriteFRB->curBlock += 1;
	recNumWrite = (recNumWrite + 1) % MaxDataBlocks;
	fprintf(stderr, "Exiting Correlator::writeToSHM_SPOTLIGHT(unsigned char *rawData, struct timeval timestamp_gps)\n");
}

void Correlator::writeToSHM(unsigned short int *rawData)
{
	fprintf(stderr, "Inside Correlator::writeToSHM(unsigned short int *rawData)\n");
	long int fetched = 0;
	ofstream warnFile;
	dataTabWrite[recNumWrite].seqnum = dataBufferWrite->cur_block;
	dataTabWrite[recNumWrite].rec = (dataBufferWrite->cur_block) % (dataBufferWrite->maxblocks);

	if (debug)
	{
		cout << endl
			 << "recNum: " << recNumWrite << endl;
		cout << "dataBufferWrite->cur_rec: " << dataBufferWrite->cur_rec << endl;
		cout << "MaxDataBuf: " << MaxDataBuf << endl;
		cout << "dataBufferWrite->cur_block: " << dataBufferWrite->cur_block << endl;
		cout << "dataBufferWrite->maxblocks: " << dataBufferWrite->maxblocks << endl;
		cout << "dataTabWrite[recNumWrite].rec: " << dataTabWrite[recNumWrite].rec << endl;
		cout << "dataTabWrite[recNumWrite].flag: " << dataTabWrite[recNumWrite].flag << endl;
		cout << "dataTabWrite[recNumWrite].seqnum: " << dataTabWrite[recNumWrite].seqnum << endl;
		cout << "dataHdrWrite->status: " << dataHdrWrite->status << endl;
		cout << "dataBufferWrite->blocksize: " << dataBufferWrite->blocksize << endl;
	}

	memcpy(dataBufferWrite->buf + dataTabWrite[recNumWrite].rec * (dataBufferWrite->blocksize) + DataOff, rawData, dataBufferWrite->blocksize - DataOff);

	dataBufferWrite->cur_rec = (recNumWrite + 1) % MaxDataBuf;
	dataBufferWrite->cur_block += 1;
	warnFile.close();
	dataTabWrite[(recNumWrite + 1) % MaxDataBuf].flag = 0;
	dataTabWrite[(recNumWrite)].flag = BufReady;
	recNumWrite = (recNumWrite + 1) % MaxDataBuf;
	fprintf(stderr, "Inside Correlator::writeToSHM(unsigned short int *rawData)\n");
}

void Correlator::writeToSHM(short int *rawData, char *header)
{
	fprintf(stderr, "Inside Correlator::writeToSHM(short int *rawData, char *header)\n");
	dataTabWrite[recNumWrite].seqnum = dataBufferWrite->cur_block;
	dataTabWrite[recNumWrite].rec = (dataBufferWrite->cur_block) % (dataBufferWrite->maxblocks);
	if (debug)
	{
		cout << endl
			 << "recNum: " << recNumWrite << endl;
		cout << "dataBufferWrite->cur_rec: " << dataBufferWrite->cur_rec << endl;
		cout << "MaxDataBuf " << MaxDataBuf << endl;
		cout << "dataBufferWrite->cur_block: " << dataBufferWrite->cur_block << endl;
		cout << "dataBufferWrite->maxblocks: " << dataBufferWrite->maxblocks << endl;
		cout << "dataTabWrite[recNumWrite].rec: " << dataTabWrite[recNumWrite].rec << endl;
		cout << "dataTabWrite[recNumWrite].flag: " << dataTabWrite[recNumWrite].flag << endl;
		cout << "dataTabWrite[recNumWrite].seqnum: " << dataTabWrite[recNumWrite].seqnum << endl;
		cout << "dataHdrWrite->status: " << dataHdrWrite->status << endl;
		cout << "dataBufferWrite->blocksize: " << dataBufferWrite->blocksize << endl;
	}
	fprintf(stderr, "Starting memcpy Buffer.\n");
	memcpy(dataBufferWrite->buf + dataTabWrite[recNumWrite].rec * (dataBufferWrite->blocksize) + DataOff, rawData, dataBufferWrite->blocksize - DataOff);
	fprintf(stderr, "Done memcpy Buffer.\n");
	fprintf(stderr, "Starting memcpy Header.\n");
	// memcpy(dataBufferWrite->buf+dataTabWrite[recNumWrite].rec*(dataBufferWrite->blocksize),header,DataOff);
	fprintf(stderr, "Done memcpy Header.\n");
	dataBufferWrite->cur_rec = (recNumWrite + 1) % MaxDataBuf;
	dataBufferWrite->cur_block += 1;
	dataTabWrite[(recNumWrite + 1) % MaxDataBuf].flag = 0;
	dataTabWrite[(recNumWrite)].flag = BufReady;
	recNumWrite = (recNumWrite + 1) % MaxDataBuf;
	fprintf(stderr, "Exiting Correlator::writeToSHM(short int *rawData, char *header)\n");
}

/*******************************************************************
 *FUNCTION: AcquireData::readFromSHM()
 *Reads from collect_psr shared memory of GSB/GWB
 *******************************************************************/
void Correlator::readFromSHM(unsigned short int *rawData)
{
	fprintf(stderr, "Inside Correlator::readFromSHM(unsigned short int *rawData)\n");
	ofstream warnFile;
	warnFile.open("realTimeWarning.gpt", ios::app);
	int flag = 0;
	while ((dataHdrRead->status == DAS_START) && (dataTabRead[recNumRead].flag & BufReady) == 0)
	{
		usleep(2000);
		if (flag == 0)
		{
			fprintf(stderr, "Waiting for DAS_START\n");
			flag = 1;
		}
	}
	if (flag == 1)
		fprintf(stderr, "Ready!\n");

	if (dataHdrRead->status != DAS_START)
	{
		if ((dataTabRead[recNumRead].flag & BufReady) == 0)
		{
			fprintf(stderr, "DAS not in START mode!!\n");
			exit(0);
		}
	}
	currentReadBlock = dataTabRead[recNumRead].seqnum;
	if (debug)
	{
		cout << endl
			 << "recNumRead: " << recNumRead << endl;
		cout << "dataBufferRead->cur_rec: " << dataBufferRead->cur_rec << endl;
		cout << "MaxDataBuf: " << MaxDataBuf << endl;
		cout << "dataBufferRead->cur_block: " << dataBufferRead->cur_block << endl;
		cout << "currentReadBlock: " << currentReadBlock << endl;
		cout << "dataBufferRead->maxblocks: " << dataBufferRead->maxblocks << endl
			 << endl;
	}

	if (dataBufferRead->cur_block - currentReadBlock >= dataBufferRead->maxblocks - 1)
	{
		warnFile << "recNum = " << recNumRead << ", Reading Sequence: " << currentReadBlock << ", Collect's Sequence: " << (dataBufferRead->cur_block - 1) << endl;
		warnFile << "Processing lagged behind..." << endl;

		fprintf(stderr, "recNumRead: %d; Reading Sequence, i.e., currentReadBlock: %ld; Collect's Sequence: %d\nRealigning...\n", recNumRead, currentReadBlock, dataBufferRead->cur_block - 1);
		recNumRead = (dataBufferRead->cur_rec - 1 - 2 + MaxDataBuf) % MaxDataBuf;
		currentReadBlock = dataTabRead[(recNumRead)].seqnum;
	}

	memcpy(rawData, dataBufferRead->buf + dataTabRead[recNumRead].rec * (dataBufferRead->blocksize) + DataOff, dataBufferRead->blocksize - DataOff);
	recNumRead = (recNumRead + 1) % MaxDataBuf;

	warnFile.close();
	fprintf(stderr, "Exiting Correlator::readFromSHM(unsigned short int *rawData)\n");
}

/*******************************************************************
 *FUNCTION: AcquireData::readFromSHM()
 *Reads from collect_psr shared memory of GSB/GWB
 *******************************************************************/
void Correlator::readFromSHM_SPOTLIGHT(unsigned char *rawData)
{
	fprintf(stderr, "Inside Correlator::readFromSHM_SPOTLIGHT(unsigned char *rawData)\n");
	ofstream warnFile;
	warnFile.open("realTimeWarning.gpt", ios::app);
	int flag = 0;
	while (currentReadBlock == dataBufferReadFRB->curRecord)
	{
		usleep(2000);
		if (flag == 0)
		{
			fprintf(stderr, "Waiting RD from FRB_SHM.\n");
			flag = 1;
		}
	}
	if (flag == 1)
		fprintf(stderr, "Ready\n");

	currentReadBlock++;
	if (debug)
	{
		cout << "dataBufferReadFRB->curRecord: " << dataBufferReadFRB->curRecord << endl;
		cout << "MaxDataBlocks: " << MaxDataBlocks << endl;
		cout << "dataBufferReadFRB->curBlock: " << dataBufferReadFRB->curBlock << endl;
		cout << "currentReadBlock: " << currentReadBlock << endl;
		cout << "recNumRead: " << recNumRead << endl;
	}

	if (dataBufferReadFRB->curBlock - currentReadBlock >= MaxDataBlocks - 1)
	{

		warnFile << "Processing lagged behind..." << endl;
		fprintf(stderr, "\nRealigning...\n");
		recNumRead = (dataBufferReadFRB->curRecord - 2 + MaxDataBlocks) % MaxDataBlocks;
		currentReadBlock = dataBufferReadFRB->curBlock;
	}
	memcpy(rawData, dataBufferReadFRB->data + ((long)DataSize * (long)recNumRead), DataSize);
	recNumRead = (recNumRead + 1) % MaxDataBlocks;

	warnFile.close();
	fprintf(stderr, "Exiting Correlator::readFromSHM_SPOTLIGHT(unsigned char *rawData)\n");
}