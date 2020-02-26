/***************************************************************************
 * Copyright 1998-2018 by authors (see AUTHORS.txt)                        *
 *                                                                         *
 *   This file is part of LuxCoreRender.                                   *
 *                                                                         *
 * Licensed under the Apache License, Version 2.0 (the "License");         *
 * you may not use this file except in compliance with the License.        *
 * You may obtain a copy of the License at                                 *
 *                                                                         *
 *     http://www.apache.org/licenses/LICENSE-2.0                          *
 *                                                                         *
 * Unless required by applicable law or agreed to in writing, software     *
 * distributed under the License is distributed on an "AS IS" BASIS,       *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
 * See the License for the specific language governing permissions and     *
 * limitations under the License.                                          *
 ***************************************************************************/

#ifndef _SLG_PMJ02_SEQUENCE_H
#define	_SLG_PMJ02_SEQUENCE_H

#include <vector>

#include <boost/thread.hpp>

#include "luxrays/core/randomgen.h"

#include "slg/slg.h"

namespace slg {

//------------------------------------------------------------------------------
// SamplePMJ
//------------------------------------------------------------------------------

class SamplePMJ {
private:
    float coordinates[2];
public:
    float& operator[](int index);
};

struct Node
{
    bool occupied, leaf;
    int min_val, max_val;
    int left, right;
};

//------------------------------------------------------------------------------
// PMJ02SampleSequenceGenerator_Pharr
//------------------------------------------------------------------------------

class PMJ02SampleSequenceGenerator_Pharr {
public:
    PMJ02SampleSequenceGenerator_Pharr(luxrays::RandomGenerator *rnd);
    void ProgressiveMultiJittered02Algorithm2D(int numberOfSamplesToGenerate, int numberOfCandidates = 10);
    void exportSampleSet(std::string &outputPath);
    std::vector<SamplePMJ> generatedSamples;

private:
    luxrays::RandomGenerator *randomNumberGenerator;
    int numberOfSamplesToGenerate;
    std::vector< std::vector<bool> > occupiedStrata;
    std::vector<int> xhalves;
    std::vector<int> yhalves;
    int numSamples;
    int numberOfCandidates;
    int grid_dim;
    std::vector<Node> x_tree, y_tree;
    std::vector<int> x_offsets, y_offsets;
    std::vector<SamplePMJ*> grid;
    int x_tree_idx, y_tree_idx;

    float min_dist;

    void extendSequenceEven(int alreadyGeneratedSamples);
    void extendSequenceOdd(int alreadyGeneratedSamples);
    void markOccupiedStrata(int alreadyGeneratedSamples);
    void markOccupiedStrata1(SamplePMJ &pt, int NN);
    void generateSamplePoint(int i, int j, int xhalf, int yhalf, float n, int N);
    bool isOccupied(SamplePMJ &pt, int NN);
    bool minDist(SamplePMJ& pt, float* min_dist);
    float generateRandomFloat();
    void initialize_x_tree( Node& node, int i, int j, int shape, int nx, int ny );
    void initialize_y_tree( Node& node, int i, int j, int shape, int nx, int ny );
    void valid_x_stratum( Node& node, int i, int j, int nx, int ny );
    void valid_y_stratum( Node& node, int i, int j, int nx, int ny );
    SamplePMJ* grid_item(int i, int j ) const { return grid[ j * grid_dim + i ]; }
    void set_grid_item( SamplePMJ* pt )
    { 
        int idx = int((*pt)[1]*grid_dim) * grid_dim + int((*pt)[0]*grid_dim);
        grid[ idx ] = pt;
    }
};

//------------------------------------------------------------------------------
// PMJ02Sequence
//------------------------------------------------------------------------------

class PMJ02Sequence {
public:
	PMJ02Sequence(luxrays::RandomGenerator *rndGen);
    PMJ02Sequence(const u_int seed);
	~PMJ02Sequence();

	void RequestSamples(const u_int size);

	float GetSample(const u_int pass, const u_int index);
	std::vector<float> GetSamples(const u_int pass);

private:
	// Generates for a single pixel index
	void RequestSamples(const u_int size, const u_int index);

	luxrays::RandomGenerator *rndGen;

	struct float2 {
		float2() {}
		float2(float xi, float yi) : x(xi), y(yi) {
		}

		float x;
		float y;
	};

	void shuffle(std::vector<SamplePMJ> points, u_int size);

	// How many samples should be generated at once
	u_int num_samples;

	// A vector with each pair of dimensions
	//	A vector with the (num_samples) 2D samples for that pixel
	std::vector<std::vector<SamplePMJ>> samplePoints;

    bool localRng;

    boost::mutex sampleGenerationMtx;

};

}

#endif	/* _SLG_PMJ02_SEQUENCE_H */