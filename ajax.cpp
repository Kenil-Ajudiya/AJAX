/*************************************************************
AJAX - The RFI Mitigation Tool for the SPOTLIGHT Project of the uGMRT

The tool has its foundation in pmon and then wpmon. It was
written from scratch keeping in mind the need for it to run
real time with the high data flow rate of the upgraded GMRT.
This was done in consultation with Yashwant Gupta. This version
uses c++ constructs such as classes to steamline the code and
to make feature expansions easier. If you're working on it
please keep in mind these two guidelines:

1> Continue using similar naming conventions.
2> This programs owes a lot of its speed to how it handles array
   access. Use pointers to access arrays instead of indexing
   wherever an array is being accessed repeatedly in a loop.
   You can find multiple examples of this throughout the code.

			-Aditya Chowdhury, 14th March 2016

************************************************************/

#include <externalLibraries.h>
#include <Information.h>
#include <BasicAnalysis.h>
#include <AcquireData.h>
#include <RFIFiltering.h>
#include <AdvancedAnalysis.h>
using namespace std;

// Benchmark variables
double timeReadData = 0.0;
double timeConvertToFloat = 0.0;
double timeBandshape = 0.0;
double timeNormalization = 0.0;
double timeZeroDM = 0.0;
double timeRFITimeStats = 0.0;
double timeRFITimeFlags = 0.0;
double timeRFITimeFlagsWrite = 0.0;
double timeRFIChanStats = 0.0;
double timeRFIChanFlag = 0.0;
double timeRFIChanFlagsWrite = 0.0;
double timeFullDMUnfilteredCalc = 0.0;
double timeFullDMCalc = 0.0;
double timeFullDMUnfilteredWrite = 0.0;
double timeFullDMWrite = 0.0;
double timeProfileCalc = 0.0;
double timeProfileUnfilteredCalc = 0.0;
double timeThread1 = 0.0;
double timeThread2 = 0.0;
double timeThread3 = 0.0;
double timeThread4 = 0.0;
double timeThread5 = 0.0;
double timeThread6 = 0.0;
double timeWaitTime = 0.0;
double fillTime;
int numberOfThreadRuns = 0;
// End benchmark variables
char *startFlags;
char *endFlags;
static bool keepRunning = true;

int nParallel = 5, nSerial = 10, nBeams = nParallel * nSerial;

class ThreadPacket
{
public:
	int noOfPol;
	int nParallel;
	int nSerial;
	BasicAnalysis **basicAnalysis;
	AcquireData *acquireData;
	BasicAnalysis **basicAnalysisWrite;
	RFIFiltering **rFIFilteringTime;
	RFIFiltering **rFIFilteringChan;
	AdvancedAnalysis **advancedAnalysis;
	AdvancedAnalysis **advancedAnalysisOld;
	ThreadPacket(int noOfPol_);
	~ThreadPacket();
	void copy(ThreadPacket *threadPacket);
	void copySelect(ThreadPacket *threadPacket);
	void freeMem();
};

ThreadPacket::ThreadPacket(int noOfPol_)
{
	fprintf(stderr, "Inside ThreadPacket() Constructor.\n");
	noOfPol = noOfPol_;
	acquireData = NULL;
	basicAnalysis = new BasicAnalysis *[noOfPol];
	basicAnalysisWrite = new BasicAnalysis *[noOfPol];
	rFIFilteringChan = new RFIFiltering *[noOfPol];
	rFIFilteringTime = new RFIFiltering *[noOfPol];
	advancedAnalysis = new AdvancedAnalysis *[noOfPol];
	advancedAnalysisOld = new AdvancedAnalysis *[noOfPol];
	for (int k = 0; k < noOfPol; k++)
	{
		basicAnalysis[k] = NULL;
		basicAnalysisWrite[k] = NULL;
		rFIFilteringTime[k] = NULL;
		rFIFilteringChan[k] = NULL;
		advancedAnalysis[k] = NULL;
		advancedAnalysisOld[k] = NULL;
	}
	fprintf(stderr, "Exiting ThreadPacket() Constructor.\n");
}

ThreadPacket::~ThreadPacket()
{
	if (acquireData != NULL)
		delete acquireData;
	delete[] basicAnalysis;
	delete[] basicAnalysisWrite;
	delete[] advancedAnalysis;
	delete[] advancedAnalysisOld;
	delete[] rFIFilteringChan;
	delete[] rFIFilteringTime;
}

void ThreadPacket::copy(ThreadPacket *threadPacket)
{
	fprintf(stderr, "Inside ThreadPacket::copy(ThreadPacket *threadPacket).\n");
	acquireData = threadPacket->acquireData;
	for (int k = 0; k < noOfPol; k++)
	{
		basicAnalysis[k] = threadPacket->basicAnalysis[k];
		rFIFilteringTime[k] = threadPacket->rFIFilteringTime[k];
		rFIFilteringChan[k] = threadPacket->rFIFilteringChan[k];
		advancedAnalysis[k] = threadPacket->advancedAnalysis[k];
	}
	fprintf(stderr, "Exiting ThreadPacket::copy(ThreadPacket *threadPacket).\n");
}

void ThreadPacket::copySelect(ThreadPacket *threadPacket)
{
	fprintf(stderr, "Inside ThreadPacket::copySelect(ThreadPacket *threadPacket).\n");
	for (int k = 0; k < noOfPol; k++)
	{
		basicAnalysisWrite[k] = threadPacket->basicAnalysis[k];
		rFIFilteringTime[k] = threadPacket->rFIFilteringTime[k];
		rFIFilteringChan[k] = threadPacket->rFIFilteringChan[k];
		advancedAnalysis[k] = threadPacket->advancedAnalysis[k];
		advancedAnalysisOld[k] = threadPacket->advancedAnalysisOld[k];
	}
	fprintf(stderr, "Exiting ThreadPacket::copySelect(ThreadPacket *threadPacket).\n");
}

void ThreadPacket::freeMem()
{
	fprintf(stderr, "Inside ThreadPacket::freeMem().\n");
	for (int i = 0; i < noOfPol; i++)
	{
		if (basicAnalysisWrite[i] != NULL)
			delete basicAnalysisWrite[i];
		if (advancedAnalysisOld[i] != NULL)
			delete advancedAnalysisOld[i];
		if (rFIFilteringChan[i] != NULL)
			delete rFIFilteringChan[i];
		if (rFIFilteringTime[i] != NULL)
			delete rFIFilteringTime[i];

		basicAnalysisWrite[i] = NULL;
		advancedAnalysisOld[i] = NULL;
		rFIFilteringChan[i] = NULL;
		rFIFilteringTime[i] = NULL;
	}
	fprintf(stderr, "Exiting ThreadPacket::freeMem().\n");
}

class Runtime
{
public:
	Information info;
	int blockIndex;
	int totalBlocks;
	char *blankTimeFlags;
	char *blankChanFlags;
	float *centralpass0, *centralpass1, *stdpass0, *stdpass1;
	float **cumulativeBandpass;
	char readDoneFlag;
	char readCompleteFlag;
	ThreadPacket **threadPacket;

	/**variables for inline mode of ajax**/
	Correlator *shmInterface;
	int nbuff;

	Runtime(Information info_);
	~Runtime();
	void initializeFiles();
	void fillPipe();
	void loopThrough();
	void closePipe();
	void quickclosePipe();
	void action(int threadPacketIndex, int actionIndex);

private:
	int nActions;
	int nParallelTemp;
	char chanFirst;
	char hasReachedEof;

	void writeAll(ThreadPacket *threadPacket);
	void testStatistics(ThreadPacket *threadPacket);
	void ioTasks(int threadPacketIndex);
	void channelTasks(int threadPacketIndex, char dosecondpass);
	void timeTasks(int threadPacketIndex, char dosecondpass);
	void floatConversionTasks(int threadPacketIndex);
	void writeFlagStats(ThreadPacket *threadPacket);
};

Runtime::Runtime(Information info_)
{
	fprintf(stderr, "Inside RunTime(Information info_) Constructor.\n");
	nActions = 4;
	nParallelTemp = nParallel;
	info = info_;
	AcquireData *acquireData = new AcquireData(info);

	info.blockSizeSamples = (acquireData->info).blockSizeSamples;
	info.blockSizeSec = (acquireData->info).blockSizeSec;

	fprintf(stderr, "Inside RunTime(info) Constructor, info.blockSizeSamples: %ld\n", info.blockSizeSamples);

	if (info.freqIntFactor > 1)
	{
		if (info.noOfChannels % info.freqIntFactor > 0)
		{
			fprintf(stderr, "\nNumber of channels should be an integer multiple of factor by which to integrate.\n");
			exit(0);
		}
		info.noOfChannels /= info.freqIntFactor;
		info.startChannel /= info.freqIntFactor;
		info.stopChannel /= info.freqIntFactor;
		info.smoothingWindowLength /= info.freqIntFactor;
		if (info.smoothingWindowLength <= 1)
			info.smoothingWindowLength = 3;
	}
	if (info.timeIntFactor > 1)
	{
		fprintf(stderr, "info.blockSizeSamples: %ld\n", info.blockSizeSamples);

		if (info.blockSizeSamples % info.timeIntFactor > 0)
		{
			fprintf(stderr, "info.blockSizeSizeSec: %f\n", info.blockSizeSec);
			info.blockSizeSamples += info.blockSizeSamples % info.timeIntFactor;
			info.blockSizeSec = info.blockSizeSamples * info.samplingInterval;
		}
		info.samplingInterval *= info.timeIntFactor;
		info.blockSizeSamples /= info.timeIntFactor;
		info.periodInSamples /= info.timeIntFactor;
	}

	shmInterface = new Correlator(acquireData->dataHdr, acquireData->dataBuffer);
	fprintf(stderr, "Initilizing write FRB_SHM_SPOTLIGHT.\n");
	shmInterface->initializeWriteSHM_SPOTLIGHT(0);
	nbuff = acquireData->nbuff;

	AcquireData::info = info;
	fprintf(stderr, "Inside RunTime(info) Constructor, info.blockSizeSamples: %ld, and acquireData->info.blockSizeSamples: %ld\n", info.blockSizeSamples, acquireData->info.blockSizeSamples);
	AcquireData::curPos = long((info.startTime / info.samplingInterval)) * info.noOfChannels * info.noOfPol * info.sampleSizeBytes;
	AcquireData::info.startTime = long(info.startTime / info.blockSizeSec) * info.blockSizeSec;

	fprintf(stderr, "Inside RunTime() constructor, info.filepath is: %s\n", info.filepath);
	info.display();

	threadPacket = new ThreadPacket *[nActions * nParallel];
	fprintf(stderr, "threadPacket = new ThreadPacket *[nActions * nParallel]; defined successfully.\n");
	fprintf(stderr, "nActions: %d; nParallel: %d.\n", nActions, nParallel);
	for (int i = 0; i < nActions * nParallel; i++)
	{
		fprintf(stderr, "Calling the ThreadPacket(info.noOfPol) Constructor.\n");
		threadPacket[i] = new ThreadPacket(info.noOfPol);
		threadPacket[i]->acquireData = acquireData; // Ask Aditya: Is this correct?? Commented out line 768 in ioTasks().
	}
	fprintf(stderr, "All threadPacket[i] initialised successfully.\n");

	centralpass0 = new float[info.noOfPol];
	centralpass1 = new float[info.noOfPol];
	stdpass0 = new float[info.noOfPol];
	stdpass1 = new float[info.noOfPol];
	cumulativeBandpass = new float *[info.noOfPol];
	for (int k = 0; k < info.noOfPol; k++)
		cumulativeBandpass[k] = new float[info.noOfChannels];
	int totalBlocksNoOff = 0;
	totalBlocks = 0;

	if (info.smoothFlagWindowLength > 0)
		info.smoothFlagWindowLength = (int)(info.smoothFlagWindowLength / info.samplingInterval);
	info.concentrationThreshold = 1.0 - info.concentrationThreshold / 100.0;

	fprintf(stderr, "Before for (int k = 0; k < info.noOfPol; k++).\n");

	for (int k = 0; k < info.noOfPol; k++) // Initializing analysis class objects.
	{
		cout << 'info.doFilteringOnly: ' << info.doFilteringOnly << endl;
		if (!info.doFilteringOnly)
			threadPacket[(nActions - 1) * nParallel]->advancedAnalysisOld[k] = new AdvancedAnalysis(info);
		threadPacket[nActions - 1]->basicAnalysis[k] = new BasicAnalysis(info);
		fprintf(stderr, "In for (int k = 0; k < info.noOfPol; k++).\n");
	}

	fprintf(stderr, "After for (int k = 0; k < info.noOfPol; k++).\n");

	blankTimeFlags = new char[info.blockSizeSamples + 1];
	blankChanFlags = new char[info.stopChannel - info.startChannel];

	for (int i = 0; i < info.blockSizeSamples + 1; i++)
		blankTimeFlags[i] = 0;

	for (int i = 0; i < info.stopChannel - info.startChannel; i++)
		blankChanFlags[i] = 0;

	blockIndex = 0;
	chanFirst = 0;
	if ((info.doTimeFlag && info.doChanFlag && (info.flagOrder == 1)) || info.doUseNormalizedData || info.doChanFlag)
		chanFirst = 1;
	fprintf(stderr, "Exiting Runtime(Information info_) Constructor\n");
}

Runtime::~Runtime()
{
	delete[] threadPacket;
	delete[] blankTimeFlags;
	delete[] blankChanFlags;
	// delete[] histogramInterval;
}

void Runtime::testStatistics(ThreadPacket *threadPacket)
{
	fprintf(stderr, "Inside Runtime::testStatistics(ThreadPacket *threadPacket).\n");
	for (int k = 0; k < info.noOfPol; k++)
	{
		double rmsPreFlag, meanPreFlag, rmsPostFlag, meanPostFlag, sum3, skewnessPostFilt;
		float *ptrZeroDM;
		float *ptrZeroDMUnfiltered;
		char *ptrFlags;
		int count;
		int l = threadPacket->basicAnalysisWrite[k]->blockLength;
		ptrZeroDM = threadPacket->basicAnalysisWrite[k]->zeroDM;
		ptrZeroDMUnfiltered = threadPacket->basicAnalysisWrite[k]->zeroDMUnfiltered;
		ptrFlags = threadPacket->rFIFilteringTime[k]->flags;
		meanPreFlag = rmsPreFlag = meanPostFlag = rmsPostFlag = sum3 = 0.0;
		count = 0;
		for (int i = 0; i < l; i++, ptrZeroDM++, ptrZeroDMUnfiltered++, ptrFlags++)
		{
			meanPreFlag += (*ptrZeroDMUnfiltered);
			rmsPreFlag += ((*ptrZeroDMUnfiltered) * (*ptrZeroDMUnfiltered));
			if (!(*ptrFlags))
			{
				meanPostFlag += (*ptrZeroDM);
				rmsPostFlag += (*ptrZeroDM) * (*ptrZeroDM);
				sum3 += (*ptrZeroDM) * (*ptrZeroDM) * (*ptrZeroDM);
				count++;
			}
		}

		meanPreFlag /= l;
		rmsPreFlag = sqrt((rmsPreFlag / l) - (meanPreFlag * meanPreFlag));

		meanPostFlag /= count;
		rmsPostFlag = sqrt((rmsPostFlag / count) - (meanPostFlag * meanPostFlag));
		skewnessPostFilt = (sum3 / count - 3 * meanPostFlag * rmsPostFlag - pow(meanPostFlag, 3)) / pow(rmsPostFlag, 3);
		ofstream statFile;
		if (info.doPolarMode)
		{
			ostringstream filename;
			filename << "stats" << k + 1 << ".gpt";
			statFile.open(filename.str().c_str(), ios::app);
		}
		else
			statFile.open("stats.gpt", ios::app);
		statFile << blockIndex - 3 << "\t" << threadPacket->rFIFilteringTime[k]->centralTendency << "\t" << meanPreFlag << "\t" << meanPostFlag << "\t";
		statFile << threadPacket->rFIFilteringTime[k]->rms << "\t" << rmsPreFlag << "\t" << rmsPostFlag << "\t";
		statFile << (threadPacket->rFIFilteringTime[k]->centralTendency) / (threadPacket->rFIFilteringTime[k]->rms) << "\t" << meanPreFlag / rmsPreFlag << "\t" << meanPostFlag / rmsPostFlag;
		statFile << "\t" << skewnessPostFilt << endl;
		statFile.close();
	}
}

void Runtime::writeFlagStats(ThreadPacket *threadPacket)
{
	fprintf(stderr, "Inside Runtime::writeFlagStats(ThreadPacket *threadPacket).\n");
	for (int k = 0; k < info.noOfPol; k++)
	{
		ofstream statFile;
		if (info.doPolarMode)
		{
			ostringstream filename;
			filename << "flag_stats" << k + 1 << ".gpt";
			statFile.open(filename.str().c_str(), ios::app);
		}
		else
			statFile.open("flag_stats.gpt", ios::app);

		char *ptrTimeFlags = threadPacket->rFIFilteringTime[k]->flags;
		char *ptrChanFlags = threadPacket->rFIFilteringChan[k]->flags;
		int l = threadPacket->basicAnalysisWrite[k]->blockLength;
		float timePercent;
		timePercent = 0;
		for (int i = 0; i < l; i++, ptrTimeFlags++)
			timePercent += *ptrTimeFlags;
		timePercent *= 100.0 / l;
		for (int i = 0; i < info.startChannel; i++)
			statFile << 100.0 << " ";
		for (int i = info.startChannel; i < info.stopChannel; i++, ptrChanFlags++)
		{
			if (!(*ptrChanFlags))
				statFile << timePercent << " ";
			else
				statFile << 100.0 << " ";
		}
		for (int i = info.stopChannel; i < info.noOfChannels; i++)
			statFile << 100.0 << " ";
		statFile << endl;
	}
}

void Runtime::initializeFiles()
{
	fprintf(stderr, "Inside Runtime::initializeFiles().\n");
	if (info.doWriteChanFlags && info.doChanFlag)
	{
		ofstream chanflagfile;
		chanflagfile.open("chanflag.gpt", ios::out | ios::trunc);
		chanflagfile.close();
	}
	if (info.doWriteTimeFlags && info.doTimeFlag)
	{
		ofstream timeflagfile;
		timeflagfile.open("timeflag.gpt", ios::out | ios::trunc);
		timeflagfile.close();
	}

	ofstream statFile;
	statFile.open("stats.gpt", ios::out | ios::trunc);
	statFile << "#window_indx\tmean_pred\tmean_pre\tmean_post\trms_pred\trms_pre\trms_post\tm/r_pred\tm/r_pre\tm/r_post" << endl;
	statFile.close();

	ofstream intensityFile;
	intensityFile.open("intensity_summary.gpt", ios::out | ios::trunc);
	intensityFile << "#First element of each line denotes the number of time samples in the block, followed by intensity in each channel" << endl;
	intensityFile.close();

	ofstream flagStatFile;
	flagStatFile.open("flag_stats.gpt", ios::out | ios::trunc);
	flagStatFile << "#Each line represents a seperate block. For the particular block, the line contains the percentage of flagged data in each channel" << endl;
	flagStatFile.close();
}

void Runtime::writeAll(ThreadPacket *threadPacket)
{
	fprintf(stderr, "Inside Runtime::writeAll(ThreadPacket *threadPacket).\n");
	unsigned char *filteredRawData = threadPacket->basicAnalysisWrite[0]->filteredRawDataChar;
	fprintf(stderr, "Start writing output Buffer; nbuff = %d\n", nbuff);
	for (int i = 0; i < nbuff; i++)
	{
		shmInterface->writeToSHM_SPOTLIGHT(filteredRawData, threadPacket->basicAnalysisWrite[0]->timestamp);
	}

	threadPacket->basicAnalysisWrite[0]->writeCurBandshape("intensity_summary.gpt");
	timeRFITimeFlagsWrite -= omp_get_wtime(); // benchmark
	if (info.doWriteTimeFlags && info.doTimeFlag)
		threadPacket->rFIFilteringTime[0]->writeFlags("timeflag.gpt");
	timeRFITimeFlagsWrite += omp_get_wtime(); // benchmark

	timeRFIChanFlagsWrite -= omp_get_wtime(); // benchmark
	if (info.doWriteChanFlags && info.doChanFlag)
		threadPacket->rFIFilteringChan[0]->writeFlags("chanflag.gpt", startFlags, info.startChannel, endFlags, info.noOfChannels - info.stopChannel);
	timeRFIChanFlagsWrite += omp_get_wtime(); // benchmark

	testStatistics(threadPacket);
	writeFlagStats(threadPacket);
}

void Runtime::fillPipe()
{
	fprintf(stderr, "Inside Runtime::fillPipe().\n");
	fillTime = omp_get_wtime(); // benchmark
	int index;
	readCompleteFlag = 0;
	readDoneFlag = 1;

	for (int k = 0; k < nActions; k++)
	{
		fprintf(stderr, "blockIndex: %d, and k: %d\n", blockIndex, k);
#pragma omp parallel for schedule(dynamic, 1)
		for (int j = 1; j < k + 1; j++)
		{
			fprintf(stderr, "j*nParallel: %d, and j: \n", j * nParallel, j);
			action(j * nParallel, j);
		}

		while (!readCompleteFlag)
			;
		if (k >= 3)
		{
			for (int ipol = 0; ipol < info.noOfPol; ipol++)
			{
				memcpy(cumulativeBandpass[ipol], threadPacket[4 * nParallel - 1]->basicAnalysis[ipol]->smoothBandshape, info.noOfChannels * sizeof(float));
			}
		}
		for (int i = 0; i < nParallel; i++)
			threadPacket[0 + i]->copySelect(threadPacket[(nActions - 1) * nParallel + i]);

		for (int i = (nActions - 1) * nParallel - 1; i >= 0; i--)
			threadPacket[i + nParallel]->copy(threadPacket[i]);

		readCompleteFlag = 0;
		readDoneFlag = 1;
		blockIndex += nParallel;
	}

	threadPacket[(nActions - 1) * nParallel]->advancedAnalysisOld = threadPacket[nParallel - 1]->advancedAnalysis;
	fillTime = (omp_get_wtime() - fillTime) / (nParallel * nActions); // benchmark
}

void Runtime::loopThrough()
{
	fprintf(stderr, "Inside Runtime::loopThrough()\n");
	double startTime, timeP, time0, time1, time2, time3, time4, time5, timeNet; // benchmark variables
	omp_set_nested(1);
	ofstream benchmarkfile;
	benchmarkfile.open("benchmark_threadtime_indv.gpt", ios::out | ios::trunc);
	while (!hasReachedEof)
	{
		numberOfThreadRuns += nParallel;
		startTime = omp_get_wtime();
		timeNet = omp_get_wtime();

#pragma omp parallel sections // filtering only mode of ajax
		{
#pragma omp section
			{
				time1 = omp_get_wtime();
				timeThread1 -= omp_get_wtime();
				action(1 * nParallel, 1);
				timeThread1 += omp_get_wtime();
				time1 = omp_get_wtime() - time1;
			}
#pragma omp section
			{
				time2 = omp_get_wtime();
				timeThread2 -= omp_get_wtime();
				action(2 * nParallel, 2);
				timeThread2 += omp_get_wtime();
				time2 = omp_get_wtime() - time2;
			}
#pragma omp section
			{
				time3 = omp_get_wtime();
				timeThread3 -= omp_get_wtime();
				action(3 * nParallel, 3);
				timeThread3 += omp_get_wtime();
				time3 = omp_get_wtime() - time3;
			}
		}

		while (!readCompleteFlag)
			;

		for (int ipol = 0; ipol < info.noOfPol; ipol++)
			memcpy(cumulativeBandpass[ipol], threadPacket[4 * nParallel - 1]->basicAnalysis[ipol]->smoothBandshape, info.noOfChannels * sizeof(float));

		for (int i = 0; i < nParallel; i++)
			threadPacket[i]->freeMem();

		for (int i = 0; i < nParallel; i++)
			threadPacket[0 + i]->copySelect(threadPacket[(nActions - 1) * nParallel + i]);

		for (int i = (nActions - 1) * nParallel - 1; i >= 0; i--)
			threadPacket[i + nParallel]->copy(threadPacket[i]);

		threadPacket[(nActions - 1) * nParallel]->advancedAnalysisOld = threadPacket[nParallel - 1]->advancedAnalysis;

		readCompleteFlag = 0;
		readDoneFlag = 1;
		if (!keepRunning)
		{
			fprintf(stderr, "\nTerminating program.\n");
			break;
		}

		if (info.doWindowDelay && info.doReadFromFile)
			while (omp_get_wtime() - startTime < info.blockSizeSec)
				;
		timeNet = omp_get_wtime() - timeNet;
		benchmarkfile << blockIndex << " " << timeP / nParallel << " " << time0 / nParallel << " " << time1 / nParallel << " " << time2 / nParallel << " " << time3 / nParallel << " " << time4 / nParallel << endl;

		blockIndex += nParallel;
	}
	benchmarkfile.close();
}

void Runtime::quickclosePipe()
{
	fprintf(stderr, "\nClosing AJAX from Runtime::quickclosePipe()\n");
	for (int k = 0; k < nParallelTemp; k++)
	{
		fprintf(stderr, "Block: %d\n", blockIndex + k - nActions * nParallel + 1);
		writeAll(threadPacket[k]);
	}
	if (info.doPolarMode)
	{
		for (int k = 0; k < info.noOfPol; k++)
		{
			// writing final bits
			ostringstream filename;
			ostringstream filenameUnfiltered;
			filename << "bandshape" << k + 1 << ".gpt";
			threadPacket[nParallelTemp - 1]->basicAnalysisWrite[k]->writeBandshape(filename.str().c_str());
			filename.str("");
			filename.clear();

			if (!info.doFilteringOnly)
			{
				filename << "profile_filtered" << k + 1 << ".gpt";
				filenameUnfiltered << "profile_unfiltered" << k + 1 << ".gpt";
				threadPacket[nParallelTemp - 1]->advancedAnalysis[k]->writeProfile(filename.str().c_str(), filenameUnfiltered.str().c_str());
			}
		}
	}
	else
	{
		threadPacket[nParallelTemp - 1]->basicAnalysisWrite[0]->writeBandshape("bandshape.gpt");
		if (!info.doFilteringOnly)
			threadPacket[nParallelTemp - 1]->advancedAnalysis[0]->writeProfile("profile_filtered.gpt", "profile_unfiltered.gpt");
	}
	for (int j = 0; j < nParallel; j++)
	{
		for (int i = 0; i < info.noOfPol; i++)
		{
			delete threadPacket[2 * nParallel + j]->basicAnalysis[i];
			delete threadPacket[3 * nParallel + j]->basicAnalysis[i];
			delete threadPacket[4 * nParallel + j]->basicAnalysis[i];

			if (!threadPacket[3 * nParallel + j]->rFIFilteringChan[i])
				delete threadPacket[3 * nParallel + j]->rFIFilteringChan[i];
			if (!threadPacket[4 * nParallel + j]->rFIFilteringChan[i])
				delete threadPacket[4 * nParallel + j]->rFIFilteringChan[i];

			if (!threadPacket[4 * nParallel + j]->rFIFilteringTime[i])
				delete threadPacket[4 * nParallel + j]->rFIFilteringTime[i];
		}
	}
	cout << endl;
}

void Runtime::closePipe()
{
	fprintf(stderr, "\nClosing AJAX from Runtime::closePipe()\n");
	for (int k = 0; k < nParallel; k++)
	{
		fprintf(stderr, "Block: %d\n", blockIndex + k - nActions * nParallel + 1);
		writeAll(threadPacket[k]);
	}
	int temp = nParallel;
	blockIndex += nParallel;
	for (int i = 1; i < nActions; i++)
	{
		nParallel = nParallelTemp;
		action(i * temp, i);
		nParallel = temp;
		for (int j = i + 1; j < nActions; j++)
		{
			action(j * nParallel, j);
		}
		if (i != nActions - 1) // The last packets is retained in memory to print the final profile from.
		{
			for (int j = 0; j < nParallel; j++)
				threadPacket[j]->freeMem();
		}
		// Packet transfers between different operations:
		for (int j = 0; j < nParallel; j++)
			threadPacket[0 + j]->copySelect(threadPacket[(nActions - 1) * nParallel + j]);

		for (int j = (nActions - 1) * nParallel - 1; j >= 0; j--)
			threadPacket[j + nParallel]->copy(threadPacket[j]);

		threadPacket[(nActions - 1) * nParallel]->advancedAnalysisOld = threadPacket[nParallel - 1]->advancedAnalysis;

		for (int k = 0; k < nParallel; k++)
		{
			fprintf(stderr, "Block: %d\n", blockIndex + k - nActions * nParallel + 1);
			writeAll(threadPacket[k]);
		}
		blockIndex += nParallel;
	}
	if (info.doPolarMode)
	{
		for (int k = 0; k < info.noOfPol; k++)
		{
			// writing final bits
			ostringstream filename;
			ostringstream filenameUnfiltered;
			filename << "bandshape" << k + 1 << ".gpt";
			threadPacket[nParallelTemp - 1]->basicAnalysisWrite[k]->writeBandshape(filename.str().c_str());
			filename.str("");
			filename.clear();

			if (!info.doFilteringOnly)
			{
				filename << "profile_filtered" << k + 1 << ".gpt";
				filenameUnfiltered << "profile_unfiltered" << k + 1 << ".gpt";
				threadPacket[nParallelTemp - 1]->advancedAnalysis[k]->writeProfile(filename.str().c_str(), filenameUnfiltered.str().c_str());
			}
		}
	}
	else
	{
		threadPacket[nParallelTemp - 1]->basicAnalysisWrite[0]->writeBandshape("bandshape.gpt");
		if (!info.doFilteringOnly)
			threadPacket[nParallelTemp - 1]->advancedAnalysis[0]->writeProfile("profile_filtered.gpt", "profile_unfiltered.gpt");
	}
	for (int j = 0; j < nParallelTemp; j++) // free'ing last packets
		threadPacket[j]->freeMem();
	cout << endl;
}

void Runtime::action(int threadPacketIndex, int actionIndex)
{
	fprintf(stderr, "Inside Runtime::action(threadPacketIndex: %d, actionIndex: %d)\n", threadPacketIndex, actionIndex);
	switch (actionIndex)
	{
	case 0:
		ioTasks(threadPacketIndex);
		break;

	case 1:
		floatConversionTasks(threadPacketIndex);
		break;

	case 2:
		channelTasks(threadPacketIndex, 0);
		timeTasks(threadPacketIndex, 0);
		channelTasks(threadPacketIndex, 1);
		break;

	case 3:
		timeTasks(threadPacketIndex, 1);
		break;
	}
}

/*******************************************************************
 *FUNCTION: ioTasks(int threadPacketIndex)
 *This function performs the tasks on thread one, i.e Reading data
 *and writing data.
 *For optimization considerations refer to document titled
 *"Multithreading Considerations for ajax"
 *******************************************************************/
void Runtime::ioTasks(int threadPacketIndex)
{
	fprintf(stderr, "Inside Runtime::ioTasks(threadPacketIndex: %d)\n", threadPacketIndex);
	double readTime, floatTime, restTime;
	char killFlag = 0;
	while (!killFlag)
	{
		ofstream benchmarkfile;
		benchmarkfile.open("benchmark_readtime.gpt", ios::app);
		// cout<<"i/o thread id:"<<sched_getcpu()<<endl;
		while (!readDoneFlag)
			;
		readDoneFlag = 0;

		for (int iBeam = 0; iBeam < nBeams; iBeam++)
		{
			ThreadPacket *thisThreadPacket = threadPacket[threadPacketIndex + (iBeam % nParallel)];
			// thisThreadPacket->acquireData = new AcquireData(blockIndex * nBeams + iBeam);

			readTime = omp_get_wtime();		 // benchmark
			timeReadData -= omp_get_wtime(); // benchmark
			thisThreadPacket->acquireData->readData(iBeam);
			readTime = omp_get_wtime() - readTime; // benchmark
			benchmarkfile << readTime << endl;

			hasReachedEof = thisThreadPacket->acquireData->hasReachedEof;
			if (hasReachedEof)
			{
				nParallelTemp = (iBeam % nParallel) + 1;
				killFlag = 1;
				break;
			}
		}
		for (int iBeam = 0; iBeam < nBeams; iBeam++)
		{
			ThreadPacket *thisThreadPacket = threadPacket[threadPacketIndex + iBeam];
			if (blockIndex * nBeams + iBeam >= nParallel * nActions)
			{
				fprintf(stderr, "Block: %d\n", blockIndex * nBeams + iBeam - nActions * nParallel + 1);
				writeAll(thisThreadPacket);
			}
		}
		readCompleteFlag = 1;
		if (!keepRunning)
			killFlag = 1;
	}
}

void Runtime::floatConversionTasks(int threadPacketIndex)
{
	fprintf(stderr, "Inside Runtime::floatConversionTasks(threadPacketIndex: %d)\n", threadPacketIndex);
	double readTime, floatTime, restTime;
	// cout<<"i/o thread id:"<<sched_getcpu()<<endl;
	ofstream benchmarkfile;
	benchmarkfile.open("benchmark_threadtime_indv.gpt", ios::app);

#pragma omp parallel for schedule(dynamic, 1)
	for (int iParallel = 0; iParallel < nParallel; iParallel++)
	{
		for (int iSerial = 0; iSerial < nSerial; iSerial++)
		{
			int i = iParallel * nSerial + iSerial;
			ThreadPacket *thisThreadPacket = threadPacket[threadPacketIndex + i];

			timeConvertToFloat -= omp_get_wtime(); // benchmark

			// thisThreadPacket->acquireData->splitRawData();

			// if (info.timeIntFactor > 1 || info.freqIntFactor > 1)
			// {
			// 	thisThreadPacket->acquireData->averageRawData(info.timeIntFactor, info.freqIntFactor);
			// }

			timeConvertToFloat += omp_get_wtime(); // benchmark
			for (int k = 0; k < info.noOfPol; k++)
			{
				thisThreadPacket->basicAnalysis[k] = new BasicAnalysis(thisThreadPacket->acquireData->splittedRawData[k], k, thisThreadPacket->acquireData->blockLength / info.timeIntFactor, thisThreadPacket->acquireData->timestamp_gps);
			}
			thisThreadPacket->basicAnalysis[0]->headerInfo = thisThreadPacket->acquireData->headerInfo;
			delete thisThreadPacket->acquireData;
			thisThreadPacket->acquireData = NULL;
		}
	}
}

/*******************************************************************
 *FUNCTION: void channelTasks(int threadPacketIndex)
 *This function performs the tasks on either thread two or three
 *depending on if channel filtering is performed before or after time
 *filtering.
 *The tasks are :Calculating bandshapes, finding outliers in bandshape
 *and normalizing data if enabled by user.
 *For optimization considerations refer to document titled
 *"Multithreading Considerations for ajax"
 *******************************************************************/
void Runtime::channelTasks(int threadPacketIndex, char dosecondpass)
{
	fprintf(stderr, "Inside Runtime::channelTasks(threadPacketIndex: %d)\n", threadPacketIndex);
#pragma omp parallel for ordered schedule(dynamic, 1)
	for (int iParallel = 0; iParallel < nParallel; iParallel++)
	{
		for (int iSerial = 0; iSerial < nSerial; iSerial++)
		{
			int t = iParallel * nSerial + iSerial;
			// cout<<"channel thread id:"<<sched_getcpu()<<endl;
			BasicAnalysis **basicAnalysis = threadPacket[threadPacketIndex + t]->basicAnalysis;
			RFIFiltering **rFIFilteringChan = threadPacket[threadPacketIndex + t]->rFIFilteringChan;
			RFIFiltering **rFIFilteringTime = threadPacket[threadPacketIndex + t]->rFIFilteringTime;
			for (int i = 0; i < info.noOfPol; i++)
			{
				if (dosecondpass == 0)
				{
					switch ((int)info.bandshapeToUse)
					{
					case 1:

						rFIFilteringChan[i] = new RFIFiltering(&(basicAnalysis[i]->bandshape[info.startChannel]), info.stopChannel - info.startChannel);
						rFIFilteringChan[i]->inputMax = basicAnalysis[i]->maxBandshape;
						rFIFilteringChan[i]->inputMin = basicAnalysis[i]->minBandshape;

						break;
					case 2:

						rFIFilteringChan[i] = new RFIFiltering(&(basicAnalysis[i]->normalizedBandshape[info.startChannel]), info.stopChannel - info.startChannel);
						rFIFilteringChan[i]->inputMax = basicAnalysis[i]->maxNormalizedBandshape;
						rFIFilteringChan[i]->inputMin = basicAnalysis[i]->minNormalizedBandshape;

						break;
					case 3:
						rFIFilteringChan[i] = new RFIFiltering(&(basicAnalysis[i]->meanToRmsBandshape[info.startChannel]), info.stopChannel - info.startChannel);
						rFIFilteringChan[i]->inputMax = basicAnalysis[i]->maxMeanToRmsBandshape;
						rFIFilteringChan[i]->inputMin = basicAnalysis[i]->minMeanToRmsBandshape;
						break;
					case 4:
						rFIFilteringChan[i] = new RFIFiltering(&(basicAnalysis[i]->meanToRmsBandshape[info.startChannel]), info.stopChannel - info.startChannel);
						rFIFilteringChan[i]->inputMax = basicAnalysis[i]->maxMeanToRmsBandshape;
						rFIFilteringChan[i]->inputMin = basicAnalysis[i]->minMeanToRmsBandshape;
						break;
					default:
						rFIFilteringChan[i] = new RFIFiltering(&(basicAnalysis[i]->bandshape[info.startChannel]), info.stopChannel - info.startChannel);
						rFIFilteringChan[i]->inputMax = basicAnalysis[i]->maxBandshape;
						rFIFilteringChan[i]->inputMin = basicAnalysis[i]->minBandshape;
						break;
					}
				}
				rFIFilteringChan[i]->generateBlankFlags();
				if ((info.doTimeFlag && info.doChanFlag && (info.flagOrder == 2 || (dosecondpass == 1))))
				{
					timeBandshape -= omp_get_wtime(); // benchmark
					basicAnalysis[i]->computeBandshape(rFIFilteringTime[i]->flags);
					basicAnalysis[i]->calculateCumulativeBandshapes();
					timeBandshape += omp_get_wtime(); // benchmark
				}
				else
				{
					timeBandshape -= omp_get_wtime(); // benchmark
					basicAnalysis[i]->computeBandshape();
					basicAnalysis[i]->calculateMeanByRMSBandshape();
					timeBandshape += omp_get_wtime(); // benchmark
				}
				timeBandshape -= omp_get_wtime(); // benchmark
#pragma omp ordered
				if (dosecondpass == 0 && info.doChanFlag)
				{
					timeRFIChanStats -= omp_get_wtime(); // benchmark
					rFIFilteringChan[i]->cutoffToRms = info.chanCutOffToRMS;
					rFIFilteringChan[i]->computeStatistics(info.chanFlagAlgo);
					timeRFIChanStats += omp_get_wtime(); // benchmark

					timeRFIChanFlag -= omp_get_wtime(); // benchmark

					rFIFilteringChan[i]->flagData();
					rFIFilteringChan[i]->generateManualFlags(info.nBadChanBlocks, info.badChanBlocks, info.startChannel);
				}
				else
				{
					timeBandshape += omp_get_wtime(); // benchmark
					if (info.doChanFlag)
					{
						timeRFIChanStats -= omp_get_wtime(); // benchmark
						rFIFilteringChan[i]->cutoffToRms = info.chanCutOffToRMS;
						rFIFilteringChan[i]->computeStatistics(info.chanFlagAlgo);
						timeRFIChanStats += omp_get_wtime(); // benchmark

						timeRFIChanFlag -= omp_get_wtime(); // benchmark

						rFIFilteringChan[i]->flagData();
						if ((int)info.bandshapeToUse == 4)
						{
							rFIFilteringChan[i]->input = &(basicAnalysis[i]->normalizedBandshape[info.startChannel]);
							rFIFilteringChan[i]->computeStatistics(info.chanFlagAlgo);
							rFIFilteringChan[i]->flagData();
						}
						rFIFilteringChan[i]->generateManualFlags(info.nBadChanBlocks, info.badChanBlocks, info.startChannel);
						timeRFIChanFlag += omp_get_wtime(); // benchmark
					}
					else
					{
						rFIFilteringChan[i]->generateBlankFlags();
						rFIFilteringChan[i]->generateManualFlags(info.nBadChanBlocks, info.badChanBlocks, info.startChannel);
					}
				}
			}
			if (info.doPolarMode && dosecondpass == 1)
			{
				// Set smooth bandshapes for cross terms
				basicAnalysis[1]->setNormalizedBandpass(basicAnalysis[0]->smoothBandshape, basicAnalysis[2]->smoothBandshape);
				basicAnalysis[3]->setNormalizedBandpass(basicAnalysis[0]->smoothBandshape, basicAnalysis[2]->smoothBandshape);
				// transfer RR OR LL flags to all
				char *ptrChanFlag1 = rFIFilteringChan[0]->flags;
				char *ptrChanFlag2 = rFIFilteringChan[2]->flags;
				for (int i = 0; i < info.stopChannel - info.startChannel; i++, ptrChanFlag1++, ptrChanFlag2++)
				{
					*ptrChanFlag1 = (*ptrChanFlag1 | *ptrChanFlag2);
				}
				for (int i = 0; i < info.noOfPol; i++)
				{
					ptrChanFlag1 = rFIFilteringChan[0]->flags;
					ptrChanFlag2 = rFIFilteringChan[i]->flags;
					for (int j = 0; j < info.stopChannel - info.startChannel; j++, ptrChanFlag1++, ptrChanFlag2++)
					{
						*ptrChanFlag2 = *ptrChanFlag1;
					}
				}
			}
		}
	}
}

/*******************************************************************
 *FUNCTION: void timeTasks(int threadPacketIndex)
 *This function performs the tasks on either thread two or three
 *depending on if time filtering is performed before or after time
 *filtering.
 *The tasks are :Calculating zeroDM & finding outliers in zeroDM.
 *For optimization considerations refer to document titled
 *"Multithreading Considerations for ajax"
 *******************************************************************/
void Runtime::timeTasks(int threadPacketIndex, char secondpass)
{
	fprintf(stderr, "Inside Runtime::timeTasks(threadPacketIndex: %d)\n", threadPacketIndex);
	float *histogramIntervalTemp = new float[info.noOfPol];
#pragma omp parallel for schedule(dynamic, 1)
	for (int iParallel = 0; iParallel < nParallel; iParallel++)
	{
		for (int iSerial = 0; iSerial < nSerial; iSerial++)
		{
			int t = iParallel * nSerial + iSerial;
			// cout<<"time thread id:"<<sched_getcpu()<<endl;
			BasicAnalysis **basicAnalysis = threadPacket[threadPacketIndex + iParallel]->basicAnalysis;
			RFIFiltering **rFIFilteringChan = threadPacket[threadPacketIndex + iParallel]->rFIFilteringChan;
			RFIFiltering **rFIFilteringTime = threadPacket[threadPacketIndex + iParallel]->rFIFilteringTime;

			for (int i = 0; i < info.noOfPol; i++)
			{
				timeNormalization -= omp_get_wtime(); // benchmark
				if (info.doUseNormalizedData && secondpass == 1)
					basicAnalysis[i]->normalizeData();
				timeNormalization += omp_get_wtime(); // benchmark
				if (secondpass == 0)
					rFIFilteringTime[i] = new RFIFiltering(basicAnalysis[i]->zeroDM, basicAnalysis[i]->blockLength);
				if (info.doChanFlag || (info.doTimeFlag && info.doChanFlag && (info.flagOrder == 1)))
				{
					timeZeroDM -= omp_get_wtime(); // benchmark
					if (blockIndex > 6 * nParallel && secondpass == 0)
					{
						basicAnalysis[i]->computeZeroDMNorm(rFIFilteringChan[i]->flags, cumulativeBandpass[i]);
					}
					else
						basicAnalysis[i]->computeZeroDM(rFIFilteringChan[i]->flags);

					timeZeroDM += omp_get_wtime(); // benchmark
				}
				else
				{
					timeZeroDM -= omp_get_wtime(); // benchmark
					basicAnalysis[i]->computeZeroDM(blankChanFlags);
					timeZeroDM += omp_get_wtime(); // benchmark
				}

				if (info.doTimeFlag)
				{
					timeRFITimeStats -= omp_get_wtime(); // benchmark
					rFIFilteringTime[i]->cutoffToRms = info.timeCutOffToRMS;
					rFIFilteringTime[i]->inputMax = basicAnalysis[i]->maxZeroDM;
					rFIFilteringTime[i]->inputMin = basicAnalysis[i]->minZeroDM;

					// rFIFilteringTime[i]->histogramInterval=histogramInterval[i];
					if (blockIndex <= 6 * nParallel)
					{

						rFIFilteringTime[i]->computeStatistics(2);
						// rFIFilteringTime[i]->histogramInterval=(basicAnalysis[i]->maxZeroDM-basicAnalysis[i]->minZeroDM)/(pow(basicAnalysis[i]->blockLength,1/3.0));
					}

					if (secondpass == 0)
					{
						// cout<<centralpass0[i]<<","<<stdpass0[i]<<endl;
						if (blockIndex > 6 * nParallel)
						{
							rFIFilteringTime[i]->centralTendency = centralpass0[i];
							rFIFilteringTime[i]->rms = stdpass0[i];
							rFIFilteringTime[i]->generateBlankFlags();
							rFIFilteringTime[i]->computeStatistics(2);
						}
						// cout<<centralpass0[i]<<","<<rFIFilteringTime[i]->centralTendency<<endl;
						// cout<<stdpass0[i]<<","<<rFIFilteringTime[i]->rms<<","<<endl;
						if (info.doMultiPointFilter)
							rFIFilteringTime[i]->multiPointFlagData(info.cutoff);
						else
							rFIFilteringTime[i]->flagData();
						int flagcnt = 0;
						for (int k = 0; k < rFIFilteringTime[i]->inputSize; k++)
						{
							flagcnt += !(rFIFilteringTime[i]->flags[k]);
						}
						rFIFilteringTime[i]->computeStatistics(info.timeFlagAlgo);
						// cout<<centralpass0[i]<<","<<rFIFilteringTime[i]->centralTendency<<endl;
						// cout<<stdpass0[i]<<","<<rFIFilteringTime[i]->rms<<","<<flagcnt<<endl;
					}
					else
					{

						if (blockIndex > 3 * nParallel)
							rFIFilteringTime[i]->rms = stdpass1[i]; // using earlier estimate of std only for determining histogram interval
						// rFIFilteringTime[i]->histogramInterval=(4.0*stdpass1[i])/(pow(basicAnalysis[i]->blockLength,1/3.0));
						rFIFilteringTime[i]->computeStatistics(info.timeFlagAlgo);
						rFIFilteringTime[i]->generateBlankFlags();
						if (info.doMultiPointFilter)
							rFIFilteringTime[i]->multiPointFlagData(info.cutoff);
						else
							rFIFilteringTime[i]->flagData();

						stdpass1[i] = rFIFilteringTime[i]->rms;

						if (info.doZeroDMSub == 1)
							basicAnalysis[i]->subtractZeroDM(rFIFilteringChan[i]->flags, rFIFilteringTime[i]->centralTendency);

						if (info.smoothFlagWindowLength > 0)
							rFIFilteringTime[i]->smoothFlags((int)info.smoothFlagWindowLength, info.concentrationThreshold);
						if (info.doFlagWidthThreshold)
							rFIFilteringTime[i]->flagWidthThreshold(info.flagWidthOn, info.flagWidthOff);
						if (info.doSubbandFiltering)
						{
							rFIFilteringTime[i]->subbandFlagging(info.nSubFilt, info.nSubDetect, info.noOfChannels, info.startChannel, info.stopChannel, basicAnalysis[i]->rawData, rFIFilteringChan[i]->flags);
						}
						if (info.doFlagWidthThreshold || info.doSubbandFiltering)
							memcpy(rFIFilteringTime[i]->flags, rFIFilteringTime[i]->flagsFilt, rFIFilteringTime[i]->inputSize);
						if (info.doFRB)
						{
							if (!info.doReplaceByMean)
								basicAnalysis[i]->getFilteredRawDataChar(rFIFilteringTime[i]->flags, rFIFilteringChan[i]->flags, 0);
							else if (info.doReplaceByMean == 1)
								basicAnalysis[i]->getFilteredRawDataChar(rFIFilteringTime[i]->flags, rFIFilteringChan[i]->flags, rFIFilteringTime[i]->centralTendency);
							else if (info.doReplaceByMean == 2)
								basicAnalysis[i]->getFilteredRawDataSmoothBshapeChar(rFIFilteringTime[i]->flags, rFIFilteringChan[i]->flags);
						}

						if (!info.doReplaceByMean && (info.doWriteFiltered2D || info.isInline))
							basicAnalysis[i]->getFilteredRawData(rFIFilteringTime[i]->flags, rFIFilteringChan[i]->flags, 0);
						if (info.doReplaceByMean == 1)
							basicAnalysis[i]->getFilteredRawData(rFIFilteringTime[i]->flags, rFIFilteringChan[i]->flags, rFIFilteringTime[i]->centralTendency);
						else if (info.doReplaceByMean == 2)
							basicAnalysis[i]->getFilteredRawDataSmoothBshape(rFIFilteringTime[i]->flags, rFIFilteringChan[i]->flags);
						timeRFITimeFlags += omp_get_wtime(); // benchmark
					}
				}
				else
				{
					rFIFilteringTime[i]->generateBlankFlags();
					if (info.doZeroDMSub == 1)
						basicAnalysis[i]->subtractZeroDM(rFIFilteringChan[i]->flags, 1);
				}
			}

			if (info.doPolarMode)
			{ // transfer RR OR LL flags to all
				char *ptrTimeFlag1 = rFIFilteringTime[0]->flags;
				char *ptrTimeFlag2 = rFIFilteringTime[2]->flags;
				for (int i = 0; i < basicAnalysis[0]->blockLength; i++, ptrTimeFlag1++, ptrTimeFlag2++)
				{
					*ptrTimeFlag1 = (*ptrTimeFlag1 | *ptrTimeFlag2);
				}
				for (int i = 0; i < info.noOfPol; i++)
				{
					ptrTimeFlag1 = rFIFilteringTime[0]->flags;
					ptrTimeFlag2 = rFIFilteringTime[i]->flags;
					for (int j = 0; j < basicAnalysis[0]->blockLength; j++, ptrTimeFlag1++, ptrTimeFlag2++)
					{
						*ptrTimeFlag2 = *ptrTimeFlag1;
					}
				}
			}
		}
	}
	RFIFiltering **rFIFilteringTime = threadPacket[threadPacketIndex + nParallel - 1]->rFIFilteringTime;
	for (int i = 0; i < info.noOfPol; i++)
	{
		for (int i = 0; i < info.noOfPol; i++)
		{
			if (secondpass == 0)
			{
				centralpass0[i] = rFIFilteringTime[i]->centralTendency;
				stdpass0[i] = rFIFilteringTime[i]->rms;
			}
			else
			{
				centralpass1[i] = rFIFilteringTime[i]->centralTendency;
				stdpass1[i] = rFIFilteringTime[i]->rms;
			}
		}
	}
	// cout<<"Pass:"<<(int)secondpass<<endl;
	// for(int i=0;i<info.noOfPol;i++)
	//{
	//	histogramInterval[i]=histogramIntervalTemp[i];
	//	cout<<histogramInterval[i]<<endl;
	// }
	// delete[] histogramIntervalTemp;
}

void intHandler(int)
{
	keepRunning = false;
}

int main(int argc, char *argv[])
{
	double totalTime = omp_get_wtime(); // benchmark
	Information info;
	info.startTime = 0, info.startBlockIndex = 0;
	info.doFilteringOnly = 1; //***** -nodedisp set as the default option. It can still be enabled by using the -Dodedisp option. *****
	info.doUseTempo2 = 0;
	info.doZeroDMSub = 0;
	info.doRunFilteredMode = 0;
	info.psrcatdbPath = NULL;
	info.isInline = 0;
	info.doFRB = 0;
	int arg = 1;
	info.meanval = 8 * 1024; // Hard-coding optimization for Band-4 FRB data

	info.doFlagWidthThreshold = 1;
	info.flagWidthOn = 3;
	info.flagWidthOff = 3;

	info.doSubbandFiltering = 1;
	info.nSubFilt = 5;
	info.nSubDetect = 3;
	info.timeIntFactor = 1;
	info.freqIntFactor = 1;

	if (argc < 2)
	{
		cout << "Not enough arguments." << endl;
		info.displayNoOptionsHelp();
		exit(0);
	}
	while ((arg = getopt(argc, argv, "f:I:s:o:m:dzgS:P:C:T:hH")) != -1)
	{
		switch (arg)
		{
		case 'T':
		{
			info.timeIntFactor = (char)info.stringToDouble(optarg);
		}
		break;

		case 'C':
		{
			info.freqIntFactor = (char)info.stringToDouble(optarg);
		}
		break;

		case 'd':
		{
			info.doFRB = 1;
		}
		break;

		case 'f':
		{
			info.filepath = optarg;
			cout << "\n\n\ninfo.filepath read from the comand line: " << info.filepath << endl;
			info.doReadFromFile = 1;
		}
		break;

		case 'o':
		{
			info.outputfilepath = optarg;
		}
		break;

		case 'S':
		{
			nSerial = info.stringToDouble(optarg);
		}
		break;

		case 'P':
		{
			nParallel = info.stringToDouble(optarg);
		}
		break;

		case 'z':
		{
			info.doZeroDMSub = 1;
		}
		break;

		case 'g':
		{
			info.doRunFilteredMode = 1;
		}
		break;

		case 'm':
		{
			info.meanval = int(info.stringToDouble(optarg));
		}
		break;

		case 'h':
		case 'H':
		default:
		{
			info.displayNoOptionsHelp();
			exit(1);
		}
		break;
		}
	}
	// if (!info.checkAjaxInputFileVersion())
	// {
	// 	cout << "Old version of ajax.in file found.\n Replacing by new version formatting.\n Old version copied to ajax.in.oldver" << endl;
	// 	cout << "Note this conversion will fail if the oldversion is not ver 1.5" << endl;
	// 	info.reformatAjaxInputFile();
	// }
	info.readAjaxInputFile();

	startFlags = new char[info.startChannel];
	endFlags = new char[info.noOfChannels - info.stopChannel];
	for (int i = 0; i < info.startChannel; i++)
		startFlags[i] = 1;
	for (int i = 0; i < info.noOfChannels - info.stopChannel; i++)
		endFlags[i] = 1;
	fprintf(stderr, "Attempting to define Runtime *runtime = new RunTime(info, nParallel).\n");
	Runtime *runtime = new Runtime(info);
	fprintf(stderr, "Runtime *runtime = new RunTime(info, nParallel) defined successfully.\n");
	// runtime->initializeFiles();

	// code to capture cltr+c termination
	struct sigaction act;
	act.sa_handler = intHandler;
	sigaction(SIGINT, &act, NULL);

	fprintf(stderr, "About to go action(0,0); \n\n");
#pragma omp parallel sections
	{
#pragma omp section
		{
			runtime->action(0, 0);
		}
#pragma omp section
		{
			runtime->fillPipe();
			double loopTime = omp_get_wtime(); // benchmark
			runtime->loopThrough();
			loopTime = omp_get_wtime() - loopTime;
			if (!keepRunning)
				runtime->quickclosePipe();
			else
				runtime->closePipe();
		}
	}

	// writing benchmark files
	int i = runtime->blockIndex;
	ofstream benchmarkfile;
	benchmarkfile.open("benchmark.gpt", ios::app);
	benchmarkfile << info.samplingInterval << ",";
	benchmarkfile << info.blockSizeSamples << ",";
	benchmarkfile << (timeReadData + timeWaitTime) / (float)(i) << ",";
	benchmarkfile << (timeConvertToFloat) / (float)(i) << ",";
	benchmarkfile << -timeWaitTime / (float)(i) << ",";
	benchmarkfile << timeBandshape / (float)(i) << ",";
	benchmarkfile << timeZeroDM / (float)(i) << ",";
	benchmarkfile << timeNormalization / (float)(i) << ",";
	benchmarkfile << timeRFITimeStats / (float)(i) << ",";
	benchmarkfile << timeRFITimeFlags / (float)(i) << ",";
	benchmarkfile << timeRFITimeFlagsWrite / (float)(i) << ",";
	benchmarkfile << timeRFIChanStats / (float)(i) << ",";
	benchmarkfile << timeRFIChanFlag / (float)(i) << ",";
	benchmarkfile << timeRFIChanFlagsWrite / (float)(i) << ",";
	benchmarkfile << timeFullDMCalc / (float)(i) << ",";
	benchmarkfile << timeFullDMWrite / (float)(i) << ",";
	benchmarkfile << timeProfileCalc / (float)(i) << ",";
	benchmarkfile << timeFullDMUnfilteredCalc / (float)(i) << ",";
	benchmarkfile << timeFullDMUnfilteredWrite / (float)(i) << ",";
	benchmarkfile << timeProfileUnfilteredCalc / (float)(i) << ",";
	benchmarkfile.close();

	totalTime = omp_get_wtime() - totalTime;

	benchmarkfile.open("benchmark_threadtime.gpt", ios::app);
	benchmarkfile << nParallel << "," << runtime->info.blockSizeSamples << "," << i << "," << timeThread1 / (float)(numberOfThreadRuns) << "," << timeThread2 / (float)(numberOfThreadRuns) << "," << timeThread3 / (float)(numberOfThreadRuns) << "," << timeThread4 / (float)(numberOfThreadRuns) << "," << timeWaitTime << "," << fillTime << "," << totalTime << endl;
	benchmarkfile.close();

	benchmarkfile.open("benchmark_fillTime.gpt", ios::app);
	benchmarkfile << nParallel << "," << fillTime << endl;
	benchmarkfile.close();
	delete runtime;
	exit(0);
}