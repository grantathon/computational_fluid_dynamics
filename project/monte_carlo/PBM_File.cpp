#include "PBM_File.hpp"

PBMFile::PBMFile() :
	imax(100),
	jmax(20),
	obstacle(0),
	fluid(1),
	fileName("")
{
	std::cout << "Default constructor of PBMFile" << std::endl;
}

PBMFile::PBMFile(int imax, int jmax, bool obstacle, bool fluid, const std::string &filename) :
	imax(imax),
	jmax(jmax),
	obstacle(obstacle),
	fluid(fluid),
	fileName(filename)
{
	std::cout << "Main constructor of PBMFile" << std::endl;
}

PBMFile::~PBMFile()
{
	std::cout << "Destructor of PBMFile" << std::endl;
}

void PBMFile::OutputStep(double height_perc, double width_perc) const
{
	if(fileName.empty())
	{
		std::cout << "Set the file name (it's empty...)" << std::endl;
		return;
	}

	// Create the file and set the header
	std::ofstream file;
	file.open(fileName.c_str());
	file << "P1" << std::endl << imax << " " << jmax << std::endl << std::endl;

	// Construct the layout of the domain
	for(int j = 1; j < jmax+1; j++)
	{
		for(int i = 1; i < imax+1; i++)
		{
			if ((j > jmax/(1/height_perc)) && (i <= imax/(1/width_perc)))
			{
				file << obstacle << " ";
			}
			else
			{
				file << fluid << " ";
			}
		}

		file << std::endl;
	}

	file.close();
}
