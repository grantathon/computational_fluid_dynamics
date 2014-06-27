#ifndef PBMFILE_HPP_
#define PBMFILE_HPP_

#include <iostream>
#include <fstream>
#include <string>

/* Class to build specific 2D domains (see "Output*" functions) for PBM file types only */
class PBMFile
{
private:
	int imax, jmax;
	bool obstacle, fluid;
	std::string fileName;

public:
	PBMFile();
	PBMFile(int imax, int jmax, bool obstacle, bool fluid, const std::string &filename);

	~PBMFile();

	// Different domain styles (add more if you like!)
	void OutputStep(double height_perc, double width_perc) const;
};

#endif /* PBMFILE_HPP_ */
