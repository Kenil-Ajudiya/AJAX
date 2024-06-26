#include "Information.h"
#include <algorithm>
using namespace std;

/*******************************************************************
 *FUNCTION: double Information::stringToDouble(const std::string& s)
 *string& s		: string to convert from
 *returns *double x 	: the converted double
 *******************************************************************/
double Information::stringToDouble(const std::string &s)
{
	istringstream i(s);
	double x;
	if (!(i >> x))
		return 0;
	return x;
}

/*******************************************************************
 *FUNCTION: void Information::readAjaxInputFile()
 *Reads data from ajax.in file.
 *******************************************************************/
void Information::readAjaxInputFile()
{
	fprintf(stderr, "\n\n\nInside Information::readAjaxInputFile()\n");
	char *fileparstr = "ajax.in";
	string line;
	int k = 1, CompleteFlag = 0;
	ifstream filepar(fileparstr, ios::in);
	if (!filepar.is_open())
	{
		fprintf(stderr, "ajax.in not found!\nA sample ajax.in file has been written to the current directory.");
		writeWpmonIn();
		exit(0);
	}

	while (getline(filepar, line))
	{
		/*******************************************************************
		 *tempStr stores the current parameter read as string
		 *This is then converted to the required data type of the parameters
		 *according to the line number.
		 *******************************************************************/
		int c = 0;
		int e = 0;
		string tempStr;
		while (line[e] != ':' && line[e] != '\0' && line[e] != '\n')
			e++;
		while (line[c] != ' ' && line[c] != '\0' && line[c] != '\t' && c < e)
		{
			tempStr += line[c];
			c++;
		}
		switch (k)
		{

		case 1:
		case 2:
		case 3:
			break;
		case 4:
		{
			modeOperation = tempStr;
			break;
		}
		case 5:
		{
			doPolarMode = char(int(stringToDouble(tempStr)));
			if (doPolarMode)
				noOfPol = 4;
			else
				noOfPol = 1;
			break;
		}
		case 6:
		{
			sampleSizeBytes = int(stringToDouble(tempStr));
			break;
		}
		case 7:
		case 8:
			break;
		case 9:
		{
			lowestFrequency = stringToDouble(tempStr);
			break;
		}
		case 10:
		{
			bandwidth = stringToDouble(tempStr);
			break;
		}
		case 11:
		{
			int temp = int(stringToDouble(tempStr));
			if (temp == -1)
				sidebandFlag = 0;
			else if (temp == 1)
				sidebandFlag = 1;
			else
				sidebandFlag = 5; // error case
			break;
		}
		case 12:
		{
			noOfChannels = int(stringToDouble(tempStr));
			break;
		}
		case 13:
		{
			samplingInterval = stringToDouble(tempStr) / 1000.0; // Converting milliseconds to seconds
			break;
		}
		case 14:
		case 15:
			break;
		case 16:
		{
			pulsarName = tempStr;
			break;
		}
		case 17:
		{
			periodInMs = stringToDouble(tempStr);
			if (periodInMs == -1)
				doFixedPeriodFolding = 0;
			else
				doFixedPeriodFolding = 1;
			break;
		}
		case 18:
		{
			dispersionMeasure = stringToDouble(tempStr);
			break;
		}
		case 19:
		case 20:
			break;
		case 21:
		{
			periodInSamples = (int)stringToDouble(tempStr);
			break;
		}
		case 22:
		{
			profileOffset = stringToDouble(tempStr);
			break;
		}
		case 23:
		{
			nCoeffPolyco = int(stringToDouble(tempStr));
			break;
		}
		case 24:
		{
			spanPolyco = int(stringToDouble(tempStr));
			break;
		}
		case 25:
		{
			maxHA = int(stringToDouble(tempStr));
			break;
		}
		case 26:
		case 27:
			break;
		case 28:
		{
			polarChanToDisplay = char(int(stringToDouble(tempStr)));
		}
		case 29:
		{
			blockSizeSec = stringToDouble(tempStr);
			break;
		}
		case 30:
		{
			doManualMode = char(int(stringToDouble(tempStr)));
			break;
		}
		case 31:
		{
			doWindowDelay = char(int(stringToDouble(tempStr)));
			break;
		}
		case 32:
		case 33:
			break;
		case 34:
		{
			startChannel = int(stringToDouble(tempStr));
			break;
		}
		case 35:
		{
			stopChannel = noOfChannels - int(stringToDouble(tempStr));
			break;
		}
		case 36:
		{
			doChanFlag = char(int(stringToDouble(tempStr)));
			break;
		}
		case 37:
		{
			bandshapeToUse = char(int(stringToDouble(tempStr))) + 1; // Plus one done to skip the now obsolete option of using mean bandshape to filter.
			break;
		}
		case 38:
		{
			chanCutOffToRMS = stringToDouble(tempStr);
			break;
		}
		case 39:
		case 40:
			break;
		case 41:
		{
			doTimeFlag = char(int(stringToDouble(tempStr)));
			break;
		}
		case 42:
		{
			doUseNormalizedData = char(int(stringToDouble(tempStr)));
			break;
		}
		case 43:
		{
			timeFlagAlgo = char(int(stringToDouble(tempStr)));
			break;
		}
		case 44:
		{
			timeCutOffToRMS = stringToDouble(tempStr);
			break;
		}
		case 45:
		case 46:
			break;
		case 47:
		{
			smoothingWindowLength = int(stringToDouble(tempStr));
			break;
		}
		case 48:
		{
			normalizationProcedure = char(int(stringToDouble(tempStr)));
			break;
		}
		case 49:
		{
			doReplaceByMean = char(int(stringToDouble(tempStr)));
			break;
		}
		case 50:
		case 51:
			break;
		case 52:
		{
			doWriteChanFlags = char(int(stringToDouble(tempStr)));
			break;
		}
		case 53:
		{
			doWriteTimeFlags = char(int(stringToDouble(tempStr)));
			break;
		}
		case 54:
		{
			doWriteFiltered2D = char(int(stringToDouble(tempStr)));
			break;
		}
		case 55:
		{
			doWriteFullDM = char(int(stringToDouble(tempStr)));
			break;
		}
		case 56:
		case 57:
			break;
		case 58:
		{
			nBadChanBlocks = char(int(stringToDouble(tempStr)));
			break;
		}
		case 59:
			break;
		case 60:
		{
			parseManFlagList(tempStr);
			break;
		}

		default:
		{
			CompleteFlag = 1;
			break;
		}
		}
		k++;
		if (CompleteFlag == 1)
			break;
	}
	// Manually filling obsolete options (removed from .in file)
	smoothFlagWindowLength = -1;
	refFrequency = 0;
	flagOrder = 1;
	chanFlagAlgo = 2;
	doMultiPointFilter = 0;
	if (flagOrder == 2 && doUseNormalizedData) // Channel filtering cannot be done after time filtering if data is normalized
	{
		fprintf(stderr, "Error in line of ajax.in.\n\
		Channel filtering cannot be done after time filtering if data is normalized.\n\
		Recommended option is to do channel filtering first.\n");
		exit(0);
	}

	errorChecks();
	fprintf(stderr, "Exiting Information::readAjaxInputFile()\n");
}

void Information::parseManFlagList(std::string &s)
{
	fprintf(stderr, "Inside Information::parseManFlagList(std::string &s)\n");
	int p = 0;
	int i = 0;
	int p1, p2, p3;
	int slen = s.length();
	badChanBlocks = new int[nBadChanBlocks * 2];
	while (p < slen)
	{
		if (i >= nBadChanBlocks)
		{
			fprintf(stderr, "Error in line 60 of ajax.in\n\
			Expected exactly %d list of bad sub-bands.\n",
					nBadChanBlocks);
			exit(0);
		}
		if (s[p] == '[')
			p1 = p;
		else
		{
			fprintf(stderr, "Error in line 60 of ajax.in\n\
			Expected [ to mark the start of a bad sub-band.\n");
			exit(0);
		}
		while (s[++p] != ',' && s[p] != '[' && s[p] != ']' && p < slen)
			;
		if (s[p] == ',')
			p2 = p;
		else
		{
			fprintf(stderr, "Error in line 60 of ajax.in\n\
			Expected a comma (,).\n");
			exit(0);
		}
		while (s[++p] != ']' && p < slen)
			;
		if (s[p] == ']')
			p3 = p;
		else
		{
			fprintf(stderr, "Error in line 60 of ajax.in\n\
			Expected ] to mark the end of a bad sub-band.\n");
			exit(0);
		}
		if (++p < slen && s[p] != ',')
		{
			fprintf(stderr, "Error in line 60 of ajax.in\n\
			Expected a comma (,) between list of sub-bands.\n");
			exit(0);
		}
		p++;
		badChanBlocks[i * 2] = int(stringToDouble(s.substr(p1 + 1, p2 - p1 - 1)));
		badChanBlocks[i * 2 + 1] = int(stringToDouble(s.substr(p2 + 1, p3 - p2 - 1)));
		i++;
	}
	if (i < nBadChanBlocks)
	{
		fprintf(stderr, "Error in line 60 of ajax.in\n\
			Expected exactly %d list of bad sub-bands.\n",
				nBadChanBlocks);

		exit(0);
	}
}

// void Information::fillParams()
// {
// 	fprintf(stderr, "Inside Information::fillParams()\n");
// 	if (doRunFilteredMode) // rerun of .gpt output file. No filtering/
// 	{
// 		doChanFlag = 0;
// 		doTimeFlag = 0;
// 		bandshapeToUse = 1;
// 		doUseNormalizedData = 0;
// 		doReplaceByMean = 0;
// 		doWriteFiltered2D = 0;
// 	}
// 	blockSizeSamples = blockSizeSec / samplingInterval;
// 	startBlockIndex = floor(startTime / blockSizeSec);
// }

void Information::errorChecks()
{
	fprintf(stderr, "Inside Information::errorChecks()\n");
	char erFlag = 0;
	if (doReadFromFile == 1)
	{
		ifstream testExistance;
		testExistance.open(filepath + to_string(0));
		fprintf(stderr, "Inside errorChecks(), filepath is: %s\n", filepath);
		if (!testExistance.is_open())
		{
			fprintf(stderr, "No file with name: %s\n", filepath);
			erFlag = 1;
		}
	}

	if (doPolarMode != 0 && doPolarMode != 1)
	{
		cout << "Error in line 5 of ajax.in:" << endl
			 << "Polarization mode must be either 0 for for Intensity or 1 for full stokes" << endl;
		erFlag = 1;
	}

	if (sampleSizeBytes != 1 && sampleSizeBytes != 2 && sampleSizeBytes != 4)
	{
		cout << "Error in line 6 of ajax.in:" << endl
			 << "The tool can only process 1 or 2 byte integer or 4 byte float." << endl;
		erFlag = 1;
	}

	if (sidebandFlag != 0 && sidebandFlag != 1)
	{
		cout << "Error in line 11 of ajax.in:" << endl
			 << "-1 for decreasing channel ordering, +1 for increasing." << endl;
		erFlag = 1;
	}

	if (profileOffset > 1 || profileOffset < 0)
	{
		cout << "Error in 22 of ajax.in:" << endl
			 << "Pulsar phase offset must be between 0 and 1" << endl;
		erFlag = 1;
	}

	if (doManualMode != 0 && doManualMode != 1)
	{
		cout << "Error in line 30 of  ajax.in:" << endl
			 << "0 for automatic update and 1 for manual update." << endl;
		erFlag = 1;
	}

	if (startChannel > noOfChannels)
	{
		cout << "Error in line 34 of  ajax.in:" << endl
			 << "Channels to skip cannot be greater than number of channels" << endl;
		erFlag = 1;
	}

	if (stopChannel < startChannel)
	{
		cout << "Error in line 35 of  ajax.in:" << endl
			 << "The sum of channels to skip in the beginning and end exceeds the number of channels" << endl;
		erFlag = 1;
	}

	if (doChanFlag != 0 && doChanFlag != 1)
	{
		cout << "Error in line 36 of  ajax.in:" << endl
			 << "0 for no frequency domain filtering and 1 for filtering on." << endl;
		erFlag = 1;
	}

	if (bandshapeToUse != 2 && bandshapeToUse != 3 && bandshapeToUse != 4)
	{
		cout << "Error in line 37 of  ajax.in:" << endl
			 << "1 for normalized bandshape; 2 for mean to rms bandshape and 3 for both." << endl;
		erFlag = 1;
	}

	if (chanCutOffToRMS <= 0)
	{
		cout << "Error in line 38 of  ajax.in:" << endl
			 << "Cutoff must be positive." << endl;
		erFlag = 1;
	}

	if (doWriteChanFlags != 0 && doWriteChanFlags != 1)
	{
		cout << "Error in line 52 of  ajax.in:" << endl
			 << "1 to write channel flags, 0 otherwise." << endl;
		erFlag = 1;
	}

	if (doTimeFlag != 0 && doTimeFlag != 1)
	{
		cout << "Error in line 41 of  ajax.in:" << endl
			 << "0 for no time domain filtering and 1 for filtering on." << endl;
		erFlag = 1;
	}

	if (doUseNormalizedData != 0 && doUseNormalizedData != 1)
	{
		cout << "Error in line 42 of  ajax.in:" << endl
			 << "1 to normalize data, 0 otherwise." << endl;
		erFlag = 1;
	}

	if (timeFlagAlgo != 2 && timeFlagAlgo != 1)
	{
		cout << "Error in line 43 of  ajax.in:" << endl
			 << "1 to use histogram based, 2 to use MAD basesd." << endl;
		erFlag = 1;
	}

	if (timeCutOffToRMS <= 0)
	{
		cout << "Error in line 44 of  ajax.in:" << endl
			 << "Cutoff must be positive." << endl;
		erFlag = 1;
	}

	if (doWriteTimeFlags != 0 && doWriteTimeFlags != 1)
	{
		cout << "Error in line 53 of  ajax.in:" << endl
			 << "1 to write time flags, 0 otherwise." << endl;
		erFlag = 1;
	}

	if ((doUseNormalizedData == 1 || (doChanFlag == 1 && bandshapeToUse == 2)) && normalizationProcedure == 1 && (smoothingWindowLength <= 0 || smoothingWindowLength >= noOfChannels))
	{
		cout << "Error in line 47 of  ajax.in:" << endl
			 << "Smooth window length must be positive and less than the number of channels." << endl;
		erFlag = 1;
	}

	if (normalizationProcedure != 2 && normalizationProcedure != 1)
	{
		cout << "Error in line 48 of  ajax.in:" << endl
			 << "1 to use cumulative smooth bandshape, 2 to use externally supplied bandshape.dat." << endl;
		erFlag = 1;
	}

	if (doReplaceByMean != 0 && doReplaceByMean != 1 && doReplaceByMean != 2)
	{
		cout << "Error in line 49 of  ajax.in:" << endl
			 << "0 to replace by zeros, 1 to replace flagged points by modal (median) value, 2 to replace by smooth bandshape" << endl;
		erFlag = 1;
	}

	if (doReplaceByMean == 1 && doUseNormalizedData != 1)
	{
		cout << "Invalid combination of choices in line 42 and 49." << endl
			 << "Data must be normalized to replace by modal (median) values.";
		erFlag = 1;
	}

	if (doReplaceByMean == 1 && doTimeFlag != 1)
	{
		cout << "Invalid combination of choices in line 41 and 49." << endl
			 << "Time filtering must be on to replace by modal (median) values.";
		erFlag = 1;
	}

	if (normalizationProcedure == 2)
	{
		ifstream testExistance;
		testExistance.open("bandshape.dat");
		if (!testExistance.is_open())
		{
			cout << "No bandshape.dat file found." << endl;
			erFlag = 1;
		}
	}

	if (doWriteFullDM != 0 && doWriteFullDM != 1)
	{
		cout << "Error in line 55 of  ajax.in:" << endl
			 << "1 to write dedispersed time series, 0 otherwise." << endl;
		erFlag = 1;
	}

	if (doWriteFiltered2D != 0 && doWriteFiltered2D != 1)
	{
		cout << "Error in line 54 of  ajax.in:" << endl
			 << "1 to write 2D time-frequency data, 0 otherwise." << endl;
		erFlag = 1;
	}

	if (nBadChanBlocks != 0)
	{
		cout << "nBadChanBlocks is: " << nBadChanBlocks << endl;
		for (int i = 0; i < nBadChanBlocks; i++)
		{
			if (badChanBlocks[i * 2] < 0 || badChanBlocks[i * 2 + 1] < 0)
			{
				cout << "Error in line 60 of ajax.in" << endl;
				cout << "In sub-band " << i + 1 << ": channel number cannot be less than zero." << endl;
				erFlag = 1;
			}
			if (badChanBlocks[i * 2] > noOfChannels || badChanBlocks[i * 2 + 1] > noOfChannels)
			{
				cout << "Error in line 60 of ajax.in" << endl;
				cout << "In sub-band " << i + 1 << ": channel number cannot be greater than the number of channels." << endl;
				erFlag = 1;
			}
			if (badChanBlocks[i * 2] >= badChanBlocks[i * 2 + 1])
			{
				cout << "Error in line 60 of ajax.in" << endl;
				cout << "In sub-band " << i + 1 << ": end channel must be strictly greater than start channel." << endl;
				erFlag = 1;
			}
		}
	}

	if (erFlag == 1)
	{
		fprintf(stderr, "Exiting Information::errorChecks() because of some error.\n");
		exit(0);
	}

	fprintf(stderr, "Exiting Information::errorChecks() because no errors found.\n");
}

void Information::calculateCutoff()
{
	fprintf(stderr, "Inside Information::CalculateCutoff()\n");
	float table[91][5] = {{1., 0.559349, 0.39004, 0.299751, 0.243506}, {1.1, 0.62054, 0.434484, 0.334715, 0.272346}, {1.2, 0.68246, 0.479754, 0.370477, 0.301926}, {1.3, 0.745061, 0.525809, 0.407005, 0.332225}, {1.4, 0.808296, 0.57261, 0.444271, 0.363221}, {1.5, 0.872124, 0.62012, 0.482245, 0.39489}, {1.6, 0.936504, 0.668302, 0.520897, 0.427211}, {1.7, 1.0014, 0.717119, 0.560199, 0.460161}, {1.8, 1.06677, 0.766538, 0.600123, 0.493716}, {1.9, 1.13259, 0.816526, 0.64064, 0.527854}, {2., 1.19882, 0.867051, 0.681723, 0.562553}, {2.1, 1.26543, 0.918084, 0.723347, 0.59779}, {2.2, 1.33241, 0.969594, 0.765486, 0.633543}, {2.3, 1.39971, 1.02156, 0.808114, 0.669793}, {2.4, 1.46733, 1.07394, 0.851208, 0.706516}, {2.5, 1.53523, 1.12673, 0.894744, 0.743694}, {2.6, 1.6034, 1.17989, 0.938702, 0.781306}, {2.7, 1.67182, 1.23341, 0.983059, 0.819334}, {2.8, 1.74046, 1.28725, 1.02779, 0.857758}, {2.9, 1.80933, 1.34142, 1.07289, 0.896561}, {3., 1.87839, 1.39587, 1.11832, 0.935725}, {3.1, 1.94763, 1.45061, 1.16408, 0.975234}, {3.2, 2.01705, 1.5056, 1.21015, 1.01507}, {3.3, 2.08662, 1.56084, 1.2565, 1.05522}, {3.4, 2.15635, 1.61631, 1.30313, 1.09567}, {3.5, 2.22621, 1.67199, 1.35002, 1.1364}, {3.6, 2.2962, 1.72787, 1.39715, 1.1774}, {3.7, 2.36631, 1.78395, 1.44451, 1.21866}, {3.8, 2.43653, 1.8402, 1.4921, 1.26016}, {3.9, 2.50685, 1.89662, 1.53989, 1.3019}, {4., 2.57727, 1.9532, 1.58788, 1.34386}, {4.1, 2.64777, 2.00992, 1.63605, 1.38602}, {4.2, 2.71836, 2.06679, 1.6844, 1.42839}, {4.3, 2.78902, 2.12378, 1.73292, 1.47094}, {4.4, 2.85975, 2.18089, 1.78159, 1.51368}, {4.5, 2.93055, 2.23812, 1.83041, 1.55659}, {4.6, 3.00141, 2.29545, 1.87936, 1.59965}, {4.7, 3.07233, 2.35289, 1.92845, 1.64288}, {4.8, 3.1433, 2.41042, 1.97767, 1.68625}, {4.9, 3.21432, 2.46804, 2.027, 1.72975}, {5., 3.28538, 2.52573, 2.07644, 1.77339}, {5.1, 3.35648, 2.58351, 2.12599, 1.81716}, {5.2, 3.42763, 2.64136, 2.17564, 1.86104}, {5.3, 3.4988, 2.69928, 2.22537, 1.90504}, {5.4, 3.57001, 2.75726, 2.2752, 1.94914}, {5.5, 3.64126, 2.8153, 2.32511, 1.99334}, {5.6, 3.71253, 2.8734, 2.3751, 2.03764}, {5.7, 3.78382, 2.93154, 2.42516, 2.08203}, {5.8, 3.85514, 2.98974, 2.47529, 2.12651}, {5.9, 3.92648, 3.04798, 2.52549, 2.17106}, {6., 3.99784, 3.10627, 2.57575, 2.2157}, {6.1, 4.06923, 3.16459, 2.62606, 2.26041}, {6.2, 4.14062, 3.22296, 2.67644, 2.30518}, {6.3, 4.21204, 3.28135, 2.72686, 2.35003}, {6.4, 4.28346, 3.33978, 2.77733, 2.39493}, {6.5, 4.3549, 3.39824, 2.82785, 2.4399}, {6.6, 4.42635, 3.45673, 2.87841, 2.48492}, {6.7, 4.49782, 3.51524, 2.92901, 2.52999}, {6.8, 4.56929, 3.57378, 2.97964, 2.57511}, {6.9, 4.64077, 3.63233, 3.03032, 2.62028}, {7., 4.71226, 3.69091, 3.08103, 2.6655}, {7.1, 4.78375, 3.74951, 3.13176, 2.71076}, {7.2, 4.85525, 3.80813, 3.18253, 2.75605}, {7.3, 4.92676, 3.86676, 3.23333, 2.80139}, {7.4, 4.99827, 3.92541, 3.28415, 2.84676}, {7.5, 5.06978, 3.98408, 3.33499, 2.89216}, {7.6, 5.1413, 4.04275, 3.38586, 2.9376}, {7.7, 5.21282, 4.10144, 3.43675, 2.98306}, {7.8, 5.28435, 4.16014, 3.48766, 3.02856}, {7.9, 5.35587, 4.21884, 3.53859, 3.07408}, {8., 5.4274, 4.27756, 3.58953, 3.11962}, {8.1, 5.49893, 4.33628, 3.64049, 3.16519}, {8.2, 5.57046, 4.39502, 3.69147, 3.21078}, {8.3, 5.64198, 4.45375, 3.74246, 3.25639}, {8.4, 5.71351, 4.5125, 3.79346, 3.30202}, {8.5, 5.78504, 4.57125, 3.84447, 3.34767}, {8.6, 5.85657, 4.63, 3.89549, 3.39334}, {8.7, 5.92809, 4.68876, 3.94653, 3.43902}, {8.8, 5.99962, 4.74752, 3.99757, 3.48472}, {8.9, 6.07114, 4.80628, 4.04862, 3.53043}, {9., 6.14266, 4.86504, 4.09968, 3.57615}, {9.1, 6.21418, 4.92381, 4.15075, 3.62189}, {9.2, 6.2857, 4.98258, 4.20182, 3.66764}, {9.3, 6.35721, 5.04135, 4.2529, 3.7134}, {9.4, 6.42872, 5.10012, 4.30398, 3.75916}, {9.5, 6.50023, 5.15889, 4.35507, 3.80494}, {9.6, 6.57174, 5.21766, 4.40616, 3.85073}, {9.7, 6.64324, 5.27643, 4.45725, 3.89652}, {9.8, 6.71474, 5.3352, 4.50835, 3.94232}, {9.9, 6.78624, 5.39397, 4.55945, 3.98813}, {10., 6.85773, 5.45273, 4.61055, 4.03394}};

	cutoff = new float[5];
	cutoff[0] = timeCutOffToRMS;

	int i1, i2;
	i1 = floor((cutoff[0] - 1) / 0.1);
	i2 = ceil((cutoff[0] - 1) / 0.1);

	for (int j = 1; j <= 4; j++)
		cutoff[j] = table[i1][j] + ((table[i2][j] - table[i1][j]) / 0.1) * (cutoff[0] - table[i1][0]);
}

/*******************************************************************
 *FUNCTION: void Information::genMJDObs()
 *Generates Modified Julian Day (MJD) at the start of observation.
 *In case of offline data it reads relevent parameters from the header
 *file and in case of online data it converts the system time to MJD.
 *******************************************************************/
void Information::genMJDObs()
{
	fprintf(stderr, "Inside Information::genMJDObs()\n");
	int YYYY, MM, DD, HH, mm, SS, ss;
	stringstream convertTime;
	long int nanoseconds;
	double sec;
	string command, linehdr;
	ostringstream headerName;
	if (doReadFromFile == 1) // Read from file
		headerName << filepath + to_string(0) << ".hdr";
	else // Real Time
		headerName << "timestamp.gpt";
	ifstream headerFile(headerName.str().c_str(), ios::in);
	if (!headerFile.is_open())
	{
		cout << headerName.str() << ":header file not found!" << endl;
		exit(1);
	}

	// Extract date and time information and check for sanity
	getline(headerFile, linehdr); // Ignore first line
	getline(headerFile, linehdr); // Get the second line and process it
	sscanf(linehdr.c_str(), "%*s %*s %d:%d:%lf", &HH, &mm, &sec);
	getline(headerFile, linehdr);
	sscanf(linehdr.c_str(), "%*s %d:%d:%d", &DD, &MM, &YYYY);
	headerFile.close();

	// Checking the sensibility of the IST Date and time Acquired:
	if (YYYY < 0 || MM < 1 || MM > 12 || DD < 1 || DD > 31 || HH < 0 || HH > 23 || mm < 0 || mm > 59 || sec < 0.0 || sec >= 60.0)
	{
		cout << "\n\nERROR ! Awkward or invalid format of IST time and Date in the raw.hdr file.\n";
		cout << "\t'raw.hdr' file should STRICTLY contain the IST time and IST date of Observation in this format\n";
		cout << "#Start time and date\n";
		cout << "IST Time: HH:mm:SS.ssss\n";
		cout << "Date: DD:MM:YYYY\n\n";
		exit(1);
	}

	if ((MM == 4 || MM == 6 || MM == 9 || MM == 11) && DD > 30)
	{
		cout << "\n\nERROR ! Given Month can not have Days > 30.\n";
		cout << "\t'raw.hdr' file should STRICTLY contain the IST time and IST date of Observation in this format\n";
		cout << "#Start time and date\n";
		cout << "IST Time: HH:mm:SS.ssss\n";
		cout << "Date: DD:MM:YYYY\n\n";
		exit(1);
	}

	// The following line is to check for whether current year is leap year
	// The logic is if(Div_by_4 And (Not_div_by_100 Or Div_by_400)) then Leap year else not.
	if ((YYYY % 4 == 0) && ((YYYY % 100 != 0) || (YYYY % 400 == 0)))
	{
		if (MM == 2 && DD > 29)
		{
			cout << "\n\nERROR ! In the file raw.hdr : Leap-year Month of Feb can not have Days > 29.\n";
			exit(1);
		}
	}
	else if (MM == 2 && DD > 28)
	{
		cout << "\n\nERROR ! In the file raw.hdr : In non-Leap-year Month of Feb can not have Days > 28.\n";
		exit(1);
	}

	// Convert IST to UTC
	convertTime.str("");
	convertTime << "TZ=\"UTC\" date \"+%d %m %Y %H %M %S %N\" -d \"" << YYYY << "-" << MM << "-" << DD << " " << HH << ":" << mm << ":" << setprecision(12) << sec << " IST\"" << "> UTC" << endl;
	system(convertTime.str().c_str());

	// Read UTC File and Get UTC Time
	ifstream UTCfile("UTC", ios::in);
	UTCfile >> DD >> MM >> YYYY >> HH >> mm >> SS >> nanoseconds;
	sec = SS + ((double)nanoseconds / 1000000000.0);
	UTCfile.close();
	system("rm -rf UTC");
	// Calculate the MJD
	stringstream calenderToMJD;
	calenderToMJD << "cal2mjd " << YYYY << " " << MM << " " << DD << " " << HH << " " << mm << " " << setprecision(12) << sec << " >calenderToMJDoutput";
	system(calenderToMJD.str().c_str());
	ifstream mjd("calenderToMJDoutput", ios::in);
	string mjdTemp;
	getline(mjd, mjdTemp);
	getline(mjd, mjdTemp);
	sscanf(mjdTemp.c_str(), "%*s %*s %lf", &MJDObs);
	system("rm -rf calenderToMJDoutput");
	mjd.close();
}

/*******************************************************************
 *FUNCTION: void Information::display()
 *Displays all run parameters to the terminal
 *******************************************************************/
void Information::display()
{
	fprintf(stderr, "Inside Information::display()\n");

	stringstream displays;
	displays << endl
			 << "ajax ver 4.6 (optimized for Band-4 FRB detection)" << endl;
	if (doFilteringOnly)
		displays << "FILTERING ONLY MODE" << endl
				 << endl;
	if (doFRB)
		displays << "INLINE MODE for FRB pipeline - read from and write to a shared memory" << endl;
	displays << "Mode of operation: " << modeOperation << endl;

	displays << endl
			 << "Lowest frequency :" << lowestFrequency << " MHz" << endl;
	displays << "Bandwidth: " << bandwidth << " MHz" << endl;

	if (sidebandFlag)
		displays << "Sideband flag: 1. Frequency of first channel is " << lowestFrequency << " MHz" << endl;
	else
		displays << "Sideband flag: -1. Frequency of first channel is " << lowestFrequency + bandwidth << " MHz" << endl;
	displays << "Number of channels: " << noOfChannels * freqIntFactor << endl;
	displays << "Sampling interval :" << samplingInterval * 1000.0 / timeIntFactor << " ms" << endl;

	displays << endl;
	if (isInline or doFRB)
		displays << "Display window size has been adjusted to accomodate integer multiples of shared memory buffer" << endl;

	displays << "Display window size :" << blockSizeSec << " sec" << endl;
	if (doManualMode && doReadFromFile)
		displays << "Window update will happen on key press." << endl;
	if (doWindowDelay && doReadFromFile)
		displays << "ajax will spend atleast " << blockSizeSec << " seconds on each window." << endl;

	displays << endl
			 << "Number of channels excluded at start: " << startChannel * freqIntFactor << endl;
	displays << "Number of channels excluded at end: " << (noOfChannels - stopChannel) * freqIntFactor << endl;

	if (nBadChanBlocks != 0)
	{
		displays << endl
				 << "Number of bad subbands to exclude: " << nBadChanBlocks << endl;
		for (int i = 0; i < nBadChanBlocks; i++)
			displays << "Sub-band " << i + 1 << " : chan # " << badChanBlocks[i * 2] << " to " << badChanBlocks[i * 2 + 1] << endl;
	}
	if (doChanFlag)
	{
		displays << endl
				 << "Channel flagging turned on. Details of procedure: " << endl;
		switch (bandshapeToUse)
		{
		case 1:
			displays << "\tMean bandshape will be used to detect RFI." << endl;
			break;
		case 2:
			displays << "\tBandhshape will be smoothened and then normalized.\n\tNormalized bandshape will be used to detect RFI." << endl;
			if (normalizationProcedure == 1)
				displays << "\t" << smoothingWindowLength << " channels will be used as window length for smoothing" << endl;
			else
				displays << "\t" << "Data from bandshape.dat will be used to normalize" << endl;
			break;
		case 3:
			displays << "\tMean to rms bandshape will be used to detect RFI." << endl;
			break;
		case 4:
			displays << "\tMean to rms bandshape as well as normalized bandshape will be used to detect RFI." << endl;
			break;
		}
		displays << "\tCutOff to RMS ratio: " << chanCutOffToRMS << endl;
	}
	else
	{
		displays << endl
				 << "No channel flagging will be performed." << endl;
	}
	if (doTimeFlag)
	{
		displays << endl
				 << "Time flagging turned on. Details of procedure: " << endl;
		if (doUseNormalizedData)
		{
			displays << "\t2-D time-frequency data will be normalized using smoothened bandshape before filtering" << endl;
			if (normalizationProcedure == 1)
				displays << "\t" << smoothingWindowLength << " channels will be used as window length for smoothing" << endl;
			else
				displays << "\t" << "Data from bandshape.dat will be used to normalize" << endl;
		}
		if (timeFlagAlgo == 1)
			displays << "\tHistogram based algorithm selected." << endl;
		else if (timeFlagAlgo == 2)
			displays << "\tMAD based algorithm selected." << endl;

		if (doMultiPointFilter)
		{
			displays << "\tMultipoint flagging ON. CutOff to rms : " << endl;
			for (int j = 0; j < 3; j++)
				cout << "\t" << cutoff[j] << " ";
			displays << endl;
		}
		else
			displays << "\tCutOff to RMS ratio: " << timeCutOffToRMS << endl;
		displays << endl
				 << "Filtered data will have a rescaled mean of " << meanval << endl;
		if (smoothFlagWindowLength > 0)
		{
			displays << "\tAny time block of length " << smoothFlagWindowLength << " seconds will be rejected if it has more than " << concentrationThreshold << " % of corrupted data" << endl;
		}
	}
	else
	{
		displays << endl
				 << "No time flagging will be performed." << endl;
	}
	if (doTimeFlag && doChanFlag)
	{
		if (doReplaceByMean)
		{
			if (timeFlagAlgo == 1)
				displays << endl
						 << "Flagged samples will be replaced by modal value of zero DM time series." << endl;
			else if (timeFlagAlgo == 2)
				displays << endl
						 << "Flagged samples will be replaced by median value of zero DM time series." << endl;
		}
		else
			displays << endl
					 << "Flagged samples will be ignored" << endl;
	}

	if (doChanFlag && doWriteChanFlags)
		displays << "Flag output to chanflag.gpt" << endl;
	if (doTimeFlag && doWriteTimeFlags)
		displays << "Flag output to timeflag.gpt" << endl;
	if (doZeroDMSub == 1)
		displays << "Appropiately scaled version of the zero DM time series will be subtracted to minimize the rms of each spectral channel." << endl;
	if (freqIntFactor > 1)
	{
		displays << endl
				 << "Channel averaging is ON." << endl;
		displays << freqIntFactor << " channels will be integrated." << endl;
		displays << "Total number of output channels is " << noOfChannels << endl;
	}

	if (timeIntFactor > 1)
	{
		displays << endl
				 << "Time averaging is ON" << endl;
		displays << timeIntFactor << " time samples will be integrated." << endl;
		displays << "Final sampling interval :" << samplingInterval * 1000.0 << " ms." << endl;
	}

	if (!doReadFromFile)
		displays << "Taking data from shared memory" << endl;
	else
		displays << "Inside display(), Raw file path: " << filepath << endl;

	cout << "displays.str().c_str()" << displays.str().c_str() << endl;

	time_t now = time(0);
	// convert now to string form
	char *dt = ctime(&now);
	ofstream f("log.gpt", ios::app);
	f << "Time of run: " << dt << endl;
	f << displays.str().c_str() << endl
	  << endl;
	f << "--------------------------------------------" << endl;
	fprintf(stderr, "Exiting Information::display()\n");
}

/*******************************************************************
 *FUNCTION: void Information::e()
 *Displays error message when no option is given while running wpmon.
 *******************************************************************/
void Information::displayNoOptionsHelp()
{
	cout << "ajax [OPTIONS]\n\n\
    [OPTIONS]\n\n\
            -f filename                 : Read from GMRT format file filename\n\
            -s start_time_in_sec        : start processing the file after skipping some time\n\
            -o output_2D_filtered_file  : path to GMRT format filtered output file\n\
            -m mean_value_of_2D_op      : mean value of output filtered file. This is the value by which the ajax normalized data is scaled by before writing out filtered file\n\
            -d {-doFRB}                 : FRB mode of ajax - make FRB_SHM, etc.\n\
            -z {-zsub}                  : optimized zero DM subtraction. Experimental feature, use with caution\n\
            -g {-gfilt}                 : turns off all filtering options - will over-ride ajax.in inputs\n\
            -S Serial                   : number of beams to process sequentially\n\
            -P Parallel                 : number of beams to process in parallel. Total number of beams = Serial * Parallel\n\
            -C chAvgFactor              : channel integration factor\n\
            -T timeAvgFactor            : time integration factor\n\
            -h, -H, {-help}, {-Help}    : print this help and exit\n";
}

/*******************************************************************
 *FUNCTION: void Information::writeWpmonIn()
 *Writes a sample ajax.in when it is not found.
 *******************************************************************/
void Information::writeWpmonIn()
{
	fprintf(stderr, "Inside Information::writeWpmonIn()\n");
	ofstream inFile("ajax.in", ios::out);
	inFile << "#*#*#ajax input file v2.0#*#*#" << endl;
	inFile << "-------------------------------------------------" << endl;
	inFile << "#****Mode of observation****#" << endl;
	inFile << "PA\t\t: Beam mode" << endl;
	inFile << "0\t\t: Polarization mode (0-> intesity, 1-> stokes data)" << endl;
	inFile << "2\t\t: Sample size of data (in bytes, usually 2)" << endl;
	inFile << "-------------------------------------------------" << endl;
	inFile << "#****Observation Paramaters****#" << endl;
	inFile << "1030\t\t: Frequency band (lowest value in Mhz)" << endl;
	inFile << "200\t\t: Bandwidth(in Mhz)" << endl;
	inFile << "-1\t\t: Sideband flag (-1-> decreasing +1-> increasing)" << endl;
	inFile << "2048\t\t: Number of channels" << endl;
	inFile << "1.31072\t\t: Sampling Interval (in ms)" << endl;
	inFile << "-------------------------------------------------" << endl;
	inFile << "#****Pulsar Parameters****#" << endl;
	inFile << "J1807-0847\t: Pulsar name" << endl;
	inFile << "-1\t\t: Pulsar period (in milliseconds)" << endl;
	inFile << "-1\t\t: DM (in pc/cc)" << endl;
	inFile << "-------------------------------------------------" << endl;
	inFile << "#****Dedispersion & Folding parameters****#" << endl;
	inFile << "-1\t\t: Number of bins in folded profile (-1 for native resolution)" << endl;
	inFile << "0\t\t: Phase offset for folding" << endl;
	inFile << "12\t\t: Number of coefficients for each polyco span (nCoeff)" << endl;
	inFile << "60\t\t: Validity of each span (nSpan in mins)" << endl;
	inFile << "12\t\t: Maximum hour angle" << endl;
	inFile << "-------------------------------------------------" << endl;
	inFile << "#****Display Parameters****#" << endl;
	inFile << "0\t\t: Polarization channel to display (0-3 or -1 for all four)" << endl;
	inFile << "1\t\t: Display window size (seconds, 0-> pulsar period)" << endl;
	inFile << "0\t\t: Update mode	    	(0-> automatic, 1-> manual)" << endl;
	inFile << "0\t\t: Time delay between window updates (0-> no delay, 1-> emulate real time)" << endl;
	inFile << "-------------------------------------------------" << endl;
	inFile << "#****Spectral line RFI mitigation options****#" << endl;
	inFile << "50\t\t: Number of channels to flag at band beginning" << endl;
	inFile << "50\t\t: Number of channels to flag at band end" << endl;
	inFile << "1\t\t: Frequency flagging options (0-> no flagging, 1-> real time calculation)" << endl;
	inFile << "1\t\t: Bandshape to use for frequency flagging (1-> normalized bandshape, 2-> mean-to-rms bandshape, 3-> Both)" << endl;
	inFile << "2.5\t\t: Threshold for frequency flagging (in units of RMS deviation)" << endl;
	inFile << "-------------------------------------------------" << endl;
	inFile << "#****Time domain impulsive RFI mitigation options****#" << endl;
	inFile << "1\t\t: Time flagging options 	(0-> no flagging, 1-> real time calculation)" << endl;
	inFile << "1\t\t: Data normalization before filtering (0-> no, 1-> yes)" << endl;
	inFile << "1\t\t: Time flagging algorithm	(1-> histogram based, 2-> MAD based)" << endl;
	inFile << "3\t\t: Threshold for time flagging (in units of RMS deviation)" << endl;
	inFile << "-------------------------------------------------" << endl;
	inFile << "#****Other options****#" << endl;
	inFile << "20\t\t: Smoothing window size for bandshape normalization (in number of channels)" << endl;
	inFile << "1\t\t: Normalization procedure (1-> cumulative smooth bandshape, 2-> externally supplied bandshape.dat)" << endl;
	inFile << "0\t\t: Replace by median values (0-> Ignore flagged samples, 1-> Replace flagged samples by window median, 2-> Replace by smooth bandshape)" << endl;
	inFile << "-------------------------------------------------" << endl;
	inFile << "#****I/O options****#" << endl;
	inFile << "0\t\t: Write channel flag file (0-> no,1-> yes)" << endl;
	inFile << "0\t\t: Write time flag file (0-> no, 1-> yes)" << endl;
	inFile << "0\t\t: write out filtered 2D raw data (0-> no, 1-> yes)" << endl;
	inFile << "0\t\t: Write out fullDM.raw	(0-> no, 1-> yes)" << endl;
	inFile << "-------------------------------------------------" << endl;
	inFile << "#****manual flagging options****#" << endl;
	inFile << "0\t\t: Number of bad channel blocks" << endl;
	inFile << "List\t\t: #in next line, example: [200,400],[1200,1400]" << endl;
	inFile << "" << endl;
}