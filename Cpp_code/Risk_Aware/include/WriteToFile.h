/*=============================================================================
 * Copyright (C) 2023 MingYi Wang
 *
 * This program is free software: you can redistribute it and/or modify it under the terms
 * of the GNU General Public License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program.
 * If not, see http://www.gnu.org/licenses/.
 *============================================================================*/

 /*==============================================================================
  * File: WriteToFile.h
  *
  * Author: MingYi Wang (based on the code by Marc Aur��le Gilles)
  *
  * Description: This file contains helper functions for writing multi-dimensional
  * Boost arrays and vectors to file
  *
  *============================================================================*/


#ifndef WRITE_TO_FILE_H
#define WRITE_TO_FILE_H

  /** ----- Libraries ----------------------------------------------------------*/
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <boost/multi_array.hpp>
#include <boost/bimap.hpp>

namespace io {

	//**This function writes the 1D Boost matrix "aMatrix" to a file with name aFilename**//

	template <class T>
	void writeToFile1D(std::string aFilename, boost::multi_array<T, 1> aMatrix) {

		// Getting dimensions of matrix
		const int aDim0 = aMatrix.shape()[0];

		// Specifying relative file path for output files
		std::string path = "output/" + aFilename;

		// Opening dataFile in write mode,
		ofstream dataFile;
		dataFile.open(path.c_str(), std::ios::app);

		if (!dataFile) {
			// if can't open the file
			std::cout << "An error occured trying to open the file" << std::endl;
			return;
		}
		else {
			// Success!
			dataFile.close(); // Closing dataFile
		}

		// Opening dataFile in write mode
		dataFile.open(path.c_str(), std::ios::binary);
		for (int i = 0; i < aDim0; ++i) {
			// Writing each element of the matrix to file
			dataFile.write(reinterpret_cast<char*> (&aMatrix[i]), sizeof(T));
		}
		dataFile.close(); // Closing dataFile
	}

	//**This function writes aMap to a file with name aFilename**//

	void writeMapToFile1D(std::string aFilename, boost::bimap<int, double> aMap) {//assuming the key is integer-type of enumeration

		typedef typename boost::bimap<int, double>::iterator iterator;
		// Specifying relative file path for output files
		std::string path = "output/" + aFilename;

		// Opening dataFile in write mode,
		ofstream dataFile;
		dataFile.open(path.c_str(), std::ios::app);

		if (!dataFile) {
			// if can't open the file
			std::cout << "An error occured trying to open the file" << std::endl;
			return;
		}
		else {
			// Success!
			dataFile.close(); // Closing dataFile
		}

		// Opening dataFile in write mode
		dataFile.open(path.c_str(), std::ios::binary);
		for (iterator iter = aMap.begin(), iend = aMap.end(); iter != iend; ++iter) {
			// Writing each element of the matrix to file
			dataFile.write(const_cast<char*>(reinterpret_cast<const char*> (&iter->right)), sizeof(double));
		}
		dataFile.close(); // Closing dataFile
	}



	//**This function writes the 2D Boost matrix "aMatrix" to a file with name aFilename**//

	template <class T>
	void writeToFile2D(std::string aFilename, boost::multi_array<T, 2> aMatrix) {

		// Getting dimensions of matrix
		const int aDim0 = aMatrix.shape()[0];
		const int aDim1 = aMatrix.shape()[1];

		// Specifying relative file path for output files
		std::string path = "output/" + aFilename;

		// Opening dataFile in write mode,
		ofstream dataFile;
		dataFile.open(path.c_str(), std::ios::app);

		if (!dataFile) {
			// if can't open the file
			std::cout << "An error occured trying to open the file" << std::endl;
			return;
		}
		else {
			// Success!
			dataFile.close(); // Closing dataFile
		}

		// Opening dataFile in write mode
		dataFile.open(path.c_str(), std::ios::binary);
		for (int i = 0; i < aDim0; ++i) {
			for (int j = 0; j < aDim1; ++j)
			{
				// Writing each element of the matrix to file
				dataFile.write(reinterpret_cast<char*> (&aMatrix[i][j]), sizeof(T));
			}
		}
		dataFile.close(); // Closing dataFile
	}


	//** This function writes the 3D Boost matrix "aMatrix" to a file with name aFilename **//

	template <class T>
	void writeToFile3D(std::string aFilename, boost::multi_array<T, 3> aMatrix) {

		// Getting dimensions of matrix
		const int aDim0 = aMatrix.shape()[0]; // budget
		const int aDim1 = aMatrix.shape()[1]; // radius
		const int aDim2 = aMatrix.shape()[2]; // theta

		// Specifying relative file path for output files
		std::string path = "output/" + aFilename;

		// Opening dataFile in write mode,
		ofstream dataFile;
		dataFile.open(path.c_str(), std::ios::app);

		if (!dataFile) {
			// if can't open the file
			std::cout << "An error occured trying to open the file" << std::endl;
			return;
		}
		else {
			// Success!
			dataFile.close();
		}

		// Opening dataFile in write mode,
		dataFile.open(path.c_str(), std::ios::binary);
		for (int i = 0; i < aDim0; ++i) {
			for (int j = 0; j < aDim1; ++j) {
				for (int k = 0; k < aDim2; ++k)
				{
					// Writing each element of the matrix to file
					dataFile.write(reinterpret_cast<char*> (&aMatrix[i][j][k]), sizeof(T));
				}
			}
		}
		dataFile.close(); // Closing file
	}



	//** This function writes the 4D Boost matrix "aMatrix" to a file with name aFilename **//

	template <class T>
	void writeToFile4D(std::string aFilename, boost::multi_array<T, 4> aMatrix) {

		// Getting dimensions of matrix
		const int aDim0 = aMatrix.shape()[0]; // mode
		const int aDim1 = aMatrix.shape()[1]; // budget
		const int aDim2 = aMatrix.shape()[2]; // radius
		const int aDim3 = aMatrix.shape()[3]; // theta

		// Specifying relative file path for output files
		std::string path = "output/" + aFilename;

		// Opening dataFile in write mode,
		ofstream dataFile;
		dataFile.open(path.c_str(), std::ios::app);

		if (!dataFile) {
			// if can't open the file
			std::cout << "An error occured trying to open the file" << std::endl;
			return;
		}
		else {
			// Success!
			dataFile.close();
		}

		// Opening dataFile in write mode,
		dataFile.open(path.c_str(), std::ios::binary);
		for (int mode = 0; mode < aDim0; ++mode) {
			for (int i = 0; i < aDim1; ++i) {
				for (int j = 0; j < aDim2; ++j) {
					for (int k = 0; k < aDim3; ++k)
					{
						// Writing each element of the matrix to file
						dataFile.write(reinterpret_cast<char*> (&aMatrix[mode][i][j][k]), sizeof(T));
					}
				}
			}
		}
		dataFile.close(); // Closing file
	}


	template <class T>
	void writeToFile4D_init(std::string aFilename, const boost::multi_array<T, 4> aMatrix, int aSlice_indx) {

		// Getting dimensions of matrix
		const int aDim0 = aMatrix.shape()[0]; // mode
		const int aDim1 = aMatrix.shape()[1]; // budget
		const int aDim2 = aMatrix.shape()[2]; // radius
		const int aDim3 = aMatrix.shape()[3]; // theta

		// Specifying relative file path for output files
		std::string path = "output/" + aFilename;

		// Opening dataFile in write mode,
		ofstream dataFile;
		dataFile.open(path.c_str(), std::ios::app);

		if (!dataFile) {
			// if can't open the file
			std::cout << "An error occured trying to open the file" << std::endl;
			return;
		}
		else {
			// Success!
			dataFile.close();
		}

		// Opening dataFile in write mode,
		dataFile.open(path.c_str(), std::ios::binary);
		for (int mode = 0; mode < aDim0; ++mode) {
			for (int j = 0; j < aDim2; ++j) {
				for (int k = 0; k < aDim3; ++k)
				{
				    short p = (short)aMatrix[mode][aSlice_indx][j][k];
				    if (k == aDim3 - 1)
				    {   //periodic in theta
				        p = (short)aMatrix[mode][aSlice_indx][j][0];
				    }
					// Writing each element of the matrix to file
					dataFile.write(reinterpret_cast<char*> (&p), sizeof(short));
				}
			}

		}
		dataFile.close(); // Closing file
	}

	template <class T>
	void writeToFile4D_init_valuefn(std::string aFilename, const boost::multi_array<T, 4> aMatrix, int aSlice_indx) {

		// Getting dimensions of matrix
		const int aDim0 = aMatrix.shape()[0]; // mode
		const int aDim1 = aMatrix.shape()[1]; // budget
		const int aDim2 = aMatrix.shape()[2]; // radius
		const int aDim3 = aMatrix.shape()[3]; // theta

		// Specifying relative file path for output files
		std::string path = "output/" + aFilename;

		// Opening dataFile in write mode,
		ofstream dataFile;
		dataFile.open(path.c_str(), std::ios::app);

		if (!dataFile) {
			// if can't open the file
			std::cout << "An error occured trying to open the file" << std::endl;
			return;
		}
		else {
			// Success!
			dataFile.close();
		}

		// Opening dataFile in write mode,
		dataFile.open(path.c_str(), std::ios::binary);
		for (int mode = 0; mode < aDim0; ++mode) {
			for (int j = 0; j < aDim2; ++j) {
				for (int k = 0; k < aDim3; ++k)
				{
					//short p = (aMatrix[mode][aSlice_indx][j][k] == -1000) ? short(aMatrix[mode][aSlice_indx][j][k]) : (short)round(aMatrix[mode][aSlice_indx][j][k] * 1e4);
					short p = (short)floor(aMatrix[mode][aSlice_indx][j][k] * 1e4);
				    if (k == aDim3 - 1)
				    {   //periodic in theta
				        p = (short)floor(aMatrix[mode][aSlice_indx][j][0] * 1e4);
				    }
					// Writing each element of the matrix to file
					dataFile.write(reinterpret_cast<char*> (&p), sizeof(short));
				}
			}

		}
		dataFile.close(); // Closing file
	}
	//** This function Appends the 4D Boost matrix "aMatrix" to a file with name aFilename **//
	template <class T>
	void AppendToFile4D(std::string aFilename, const boost::multi_array<T, 4> aMatrix, int aSlice_indx) {

		// Getting dimensions of matrix
		const int aDim0 = aMatrix.shape()[0]; // mode
		const int aDim1 = aMatrix.shape()[1]; // budget
		const int aDim2 = aMatrix.shape()[2]; // radius
		const int aDim3 = aMatrix.shape()[3]; // theta

		// Specifying relative file path for output files
		std::string path = "output/" + aFilename;

		// Opening dataFile in write mode,
		ofstream dataFile;
		dataFile.open(path.c_str(), std::ios::binary | std::ios::app);
		for (int mode = 0; mode < aDim0; ++mode) {
			for (int j = 0; j < aDim2; ++j) {
				for (int k = 0; k < aDim3; ++k)
				{
					short p = (short)aMatrix[mode][aSlice_indx][j][k];
				    if (k == aDim3 - 1)
				    {   //periodic in theta
				        p = (short)aMatrix[mode][aSlice_indx][j][0];
				    }
					// Writing each element of the matrix to file
					dataFile.write(reinterpret_cast<char*> (&p), sizeof(short));
				}
			}
		}
		dataFile.close(); // Closing file
	}

	//** This function Appends the 4D Boost matrix "aMatrix" to a file with name aFilename **//
	template <class T>
	void AppendToFile4D_valuefn(std::string aFilename, const boost::multi_array<T, 4> aMatrix, int aSlice_indx) {

		// Getting dimensions of matrix
		const int aDim0 = aMatrix.shape()[0]; // mode
		const int aDim1 = aMatrix.shape()[1]; // budget
		const int aDim2 = aMatrix.shape()[2]; // radius
		const int aDim3 = aMatrix.shape()[3]; // theta

		// Specifying relative file path for output files
		std::string path = "output/" + aFilename;

		// Opening dataFile in write mode,
		ofstream dataFile;
		dataFile.open(path.c_str(), std::ios::binary | std::ios::app);
		for (int mode = 0; mode < aDim0; ++mode) {
			for (int j = 0; j < aDim2; ++j) {
				for (int k = 0; k < aDim3; ++k)
				{
					//short p = (aMatrix[mode][aSlice_indx][j][k] == -1000) ? short(aMatrix[mode][aSlice_indx][j][k]) : (short)round(aMatrix[mode][aSlice_indx][j][k] * 1e4);
					short p = (short)floor(aMatrix[mode][aSlice_indx][j][k] * 1e4);
				    if (k == aDim3 - 1)
				    {   //periodic in theta
				        p = (short)floor(aMatrix[mode][aSlice_indx][j][0] * 1e4);
				    }
					// Writing each element of the matrix to file
					dataFile.write(reinterpret_cast<char*> (&p), sizeof(short));
				}
			}
		}
		dataFile.close(); // Closing file
	}

	//** This function writes a 1D vector "aVec" to a file with name aFilename **//

	template <class T>
	void writeVectorToFile(std::string aFilename, std::vector<T> aVec) {

		// Specifying relative file path for output files
		std::string path = "output/" + aFilename;

		// Opening dataFile in write mode,
		ofstream dataFile;
		dataFile.open(path.c_str(), std::ios::app);

		if (!dataFile) {
			// if can't open the file
			std::cout << "An error occured trying to open the file" << std::endl;
			return;
		}
		else {
			// Success!
			dataFile.close();
		}

		// Opening dataFile in write mode
		dataFile.open(path.c_str(), std::ios::binary);
		for (int i = 0; i < aVec.size(); ++i) {
			// Write each element of vector to file
			dataFile.write(reinterpret_cast<char*> (&aVec[i]), sizeof(T));
		}
		dataFile.close(); // Closing file
	}


	template <class T>
	void writeToFile4D_init_out(std::string aFilename, const boost::multi_array<T, 4> aMatrix, int aSlice_indx) {

		// Getting dimensions of matrix
		const int aDim0 = aMatrix.shape()[0]; // mode
		const int aDim1 = aMatrix.shape()[1]; // budget
		const int aDim2 = aMatrix.shape()[2]; // radius
		const int aDim3 = aMatrix.shape()[3]; // theta

		// Specifying relative file path for output files
		std::string path = "output/" + aFilename;

		// Opening dataFile in write mode,
		ofstream dataFile;
		dataFile.open(path.c_str(), std::ios::app);

		if (!dataFile) {
			// if can't open the file
			std::cout << "An error occured trying to open the file" << std::endl;
			return;
		}
		else {
			// Success!
			dataFile.close();
		}

		// Opening dataFile in write mode,
		dataFile.open(path.c_str(), std::ios::binary);
		for (int mode = 0; mode < aDim0; ++mode) {
			for (int j = 0; j < aDim2; ++j) {
				for (int k = 0; k < aDim3; ++k)
				{
					//short p = (short)aMatrix[mode][aSlice_indx][j][k];
					//if (k == aDim3 - 1)
					//{   //periodic in theta
					//	p = (short)aMatrix[mode][aSlice_indx][j][0];
					//}
					// Writing each element of the matrix to file
					//dataFile.write(reinterpret_cast<char*> (&p), sizeof(short));

					if (k == aDim3 - 1)
					{   //periodic in theta
						dataFile.write((char*) (&aMatrix[mode][aSlice_indx][j][0]), sizeof(T));
					}
					else {
						// Writing each element of the matrix to file
						dataFile.write((char*) (&aMatrix[mode][aSlice_indx][j][k]), sizeof(T));
					}

				}
			}
		}
		dataFile.close(); // Closing file
	}

	template <class T>
	void AppendToFile4D_out(std::string aFilename, const boost::multi_array<T, 4> aMatrix, int aSlice_indx) {

		// Getting dimensions of matrix
		const int aDim0 = aMatrix.shape()[0]; // mode
		const int aDim1 = aMatrix.shape()[1]; // budget
		const int aDim2 = aMatrix.shape()[2]; // radius
		const int aDim3 = aMatrix.shape()[3]; // theta

		// Specifying relative file path for output files
		std::string path = "output/" + aFilename;

		// Opening dataFile in write mode,
		ofstream dataFile;
		dataFile.open(path.c_str(), std::ios::binary | std::ios::app);
		for (int mode = 0; mode < aDim0; ++mode) {
			for (int j = 0; j < aDim2; ++j) {
				for (int k = 0; k < aDim3; ++k)
				{
					if (k == aDim3 - 1)
					{   //periodic in theta
						dataFile.write((char*)(&aMatrix[mode][aSlice_indx][j][0]), sizeof(T));
					}
					else {
						// Writing each element of the matrix to file
						dataFile.write((char*) (&aMatrix[mode][aSlice_indx][j][k]), sizeof(T));
					}
				}
			}
		}
		dataFile.close(); // Closing file
	}

}

#endif // !WRITE_TO_FILE_H
