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

#include <math.h>

#include "slg/samplers/pmj02sequence.h"

using namespace std;
using namespace luxrays;
using namespace slg;

//------------------------------------------------------------------------------
// SamplePMJ
//------------------------------------------------------------------------------

float& SamplePMJ::operator[](int index)
{
    return coordinates[index];
}

bool isPowerOfFour(unsigned int n) 
{ 
    return n !=0 && ((n&(n-1)) == 0) && !(n & 0xAAAAAAAA); 
} 

PMJ02SampleSequenceGenerator_Pharr::PMJ02SampleSequenceGenerator_Pharr(luxrays::RandomGenerator *rnd): randomNumberGenerator(rnd) {
}

void PMJ02SampleSequenceGenerator_Pharr::ProgressiveMultiJittered02Algorithm2D(int numberOfSamplesToGenerate,  
                                                                               int numberOfCandidates) {
    this->numberOfSamplesToGenerate = numberOfSamplesToGenerate;
    int arraySize = 4*numberOfSamplesToGenerate;
    generatedSamples.resize(arraySize);
    int next_pow_4_samples = numberOfSamplesToGenerate;
    if ( !isPowerOfFour( next_pow_4_samples ) )
    {
        next_pow_4_samples = numberOfSamplesToGenerate >> 1 << 2;
        if ( !isPowerOfFour( next_pow_4_samples ) )
            next_pow_4_samples = next_pow_4_samples << 1;
    }

    int bits = int(log2(next_pow_4_samples << 1 ));

    occupiedStrata.reserve(bits);
    for(int i = 0; i < bits;i++)
    {
       std::vector< bool > temp( next_pow_4_samples, false );
       occupiedStrata.push_back( temp );
    }
    xhalves.resize(4 * next_pow_4_samples);
    yhalves.resize(4*next_pow_4_samples);

    grid_dim = int(ceil(sqrt(next_pow_4_samples)));
    std::vector<SamplePMJ*> grid_temp(next_pow_4_samples, nullptr);
    grid.swap( grid_temp );

    min_dist = next_pow_4_samples;

    generatedSamples[0][0] = generateRandomFloat();
    generatedSamples[0][1] = generateRandomFloat();
    set_grid_item( &(generatedSamples[0]) );

    {
       std::vector< Node > x_temp( pow( 2, bits ), Node() );
       x_tree.swap( x_temp );
       std::vector< Node > y_temp( pow( 2, bits ), Node() );
       y_tree.swap( y_temp );
    }
    
    {
       std::vector< int > x_temp( pow( 2, bits - 1 ), -1 );
       x_offsets.swap( x_temp );
       std::vector< int > y_temp( pow( 2, bits - 1 ), -1);
       y_offsets.swap( y_temp );
    }
    
    numSamples = 1;
    this->numberOfCandidates = numberOfCandidates;
    int currentlyGeneratedSamples = 1;
    while(currentlyGeneratedSamples < numberOfSamplesToGenerate)
    {
        extendSequenceEven(currentlyGeneratedSamples);
        extendSequenceOdd(currentlyGeneratedSamples*2);
        currentlyGeneratedSamples*=4;
    }
}

void PMJ02SampleSequenceGenerator_Pharr::exportSampleSet(std::string &outputPath) {
    std::ofstream outputStream(outputPath);
    outputStream.precision(16);
    if(outputStream) {
        for (int i = 0; i < numberOfSamplesToGenerate; i++) {
            outputStream << generatedSamples[i][0];
            outputStream << " ";
            outputStream << generatedSamples[i][1];
            outputStream << "\n";
        }
        outputStream.close();
    }
    else
        std::cout << "ERREUR d'ouverture de fichier" << std::endl;
}

void PMJ02SampleSequenceGenerator_Pharr::extendSequenceEven(int alreadyGeneratedSamples) {

    float n = sqrt((float)alreadyGeneratedSamples);
    markOccupiedStrata(alreadyGeneratedSamples);
    for(int s = 0; s < alreadyGeneratedSamples; s++)
    {
        SamplePMJ oldSample = generatedSamples[s];
        int i = (int)(n*oldSample[0]);
        int j = (int)(n*oldSample[1]);
        unsigned int xhalf = (int)(2.0*(n*oldSample[0] - i));
        unsigned int yhalf = (int)(2.0*(n*oldSample[1] - j));
        xhalf = 1-xhalf;
        yhalf = 1-yhalf;
        generateSamplePoint(i,j,xhalf,yhalf,n,alreadyGeneratedSamples);
    }
}

void PMJ02SampleSequenceGenerator_Pharr::extendSequenceOdd(int alreadyGeneratedSamples) {

    float n = sqrt((float)alreadyGeneratedSamples/2.0f);
    markOccupiedStrata(alreadyGeneratedSamples);
    for(int s  = 0; s < alreadyGeneratedSamples/2; s++)
    {
        SamplePMJ oldSample = generatedSamples[s];
        int i = (int)(n*oldSample[0]);
        int j = (int)(n*oldSample[1]);
        unsigned int xhalf = (int)(2.0*(n*oldSample[0] - i));
        unsigned int yhalf = (int)(2.0*(n*oldSample[1] - j));
        if(generateRandomFloat() > 0.5)
            xhalf = 1-xhalf;
        else
            yhalf = 1-yhalf;
        xhalves[s] = xhalf;
        yhalves[s] = yhalf;
        generateSamplePoint(i,j,xhalf,yhalf,n,alreadyGeneratedSamples);
    }
    for(int s = 0; s < alreadyGeneratedSamples/2; s++)
    {
        SamplePMJ oldSample = generatedSamples[s];
        int i = (int)(n*oldSample[0]);
        int j = (int)(n*oldSample[1]);
        int xhalf = 1-xhalves[s];
        int yhalf = 1-yhalves[s];
        generateSamplePoint(i,j,xhalf,yhalf,n,alreadyGeneratedSamples);
    }

}

void PMJ02SampleSequenceGenerator_Pharr::markOccupiedStrata(int alreadyGeneratedSamples) {
    int NN = 2*alreadyGeneratedSamples;
    for(int i = 0; i <= log2(NN); i++)
        for(int j = 0; j < NN; j++)
            occupiedStrata[i][j] = false;

    for(int s = 0; s < alreadyGeneratedSamples; s++)
    {
        markOccupiedStrata1(generatedSamples[s], NN);
    }
}

void PMJ02SampleSequenceGenerator_Pharr::markOccupiedStrata1(SamplePMJ &pt, int NN) {
    int shape = 0;
    int xdivs = NN;
    int ydivs = 1;
    do{
        int xstratum = (int)(xdivs * pt[0]);
        int ystratum = (int)(ydivs * pt[1]);
        occupiedStrata[shape][ystratum*xdivs+xstratum] = true;
        shape ++;
        xdivs /= 2;
        ydivs *= 2;
    }while(xdivs != 0);
}

void PMJ02SampleSequenceGenerator_Pharr::generateSamplePoint(int i, int j, int xhalf, int yhalf, float n, int N) {

    int NN = 2*N;
    int lg2NN = int(log2(NN));

    int nx = int(pow(2, ceil(lg2NN/2.)) );
    int nx_2 = nx/n;
    bool is_even = lg2NN % 2 == 0;
    int ny = is_even ? nx : int(pow(2, floor(lg2NN/2.)) );
    int ny_2 = ny/n;


    float bestDist = 0.0;
    SamplePMJ pt = SamplePMJ();
    x_tree.clear();
    y_tree.clear();
    x_offsets.clear();
    y_offsets.clear();

    Node& node_x = x_tree[0];
    x_tree_idx = 0;
    initialize_x_tree( node_x, i*nx_2 + xhalf, j*ny_2 + (is_even ? yhalf : 0), 
                       floor(lg2NN/2.), nx, ny );
    valid_x_stratum( node_x, i*nx_2 + xhalf, j*ny_2 + (is_even ? yhalf : 0), nx, ny );

    Node& node_y = y_tree[0];
    y_tree_idx = 0;
   initialize_y_tree( node_y, i*ny_2 + (is_even ? xhalf : 0), j*nx_2 + yhalf,
                      ceil( lg2NN / 2. ), ny, nx );
   valid_y_stratum( node_y, i*ny_2 + (is_even ? xhalf : 0), j*nx_2 + yhalf, ny, nx );

    SamplePMJ candpt = SamplePMJ();

    int i_pt = x_offsets[ 0 ];
    int j_pt = y_offsets[ 0 ];

      candpt[ 0 ] = ( i_pt + generateRandomFloat() ) / NN;
      candpt[ 1 ] = ( j_pt + generateRandomFloat() ) / NN;

    int t = 1;
    do
    {

        float d;
        bool has_neighbors = minDist(candpt, &d );
        if ( has_neighbors )
        {
            if(d>bestDist)
            {
                bestDist = d;
                pt = candpt;
            }
        }
        else
            pt = candpt;

        if ( ++t < numberOfCandidates )
        {
            candpt[0] = (i_pt+generateRandomFloat())/NN;
            candpt[1] = (j_pt+generateRandomFloat())/NN;
        }
    } while ( t < numberOfCandidates );

    if ( bestDist > 0 && bestDist < min_dist )
        min_dist = bestDist;

    markOccupiedStrata1(pt, NN);
    generatedSamples[numSamples] = pt;
    set_grid_item( &(generatedSamples[numSamples]) );

    numSamples++;
}

bool PMJ02SampleSequenceGenerator_Pharr::isOccupied(SamplePMJ &pt, int NN) {
    int shape = 0;
    int xdivs = NN;
    int ydivs = 1;
    do{
        int xstratum = (int)(xdivs * pt[0]);
        int ystratum = (int)(ydivs * pt[1]);
        if(occupiedStrata[shape][ystratum*xdivs+xstratum])
            return true;
        shape++;
        xdivs/=2;
        ydivs*=2;
    }while(xdivs != 0);
    return false;
}

bool PMJ02SampleSequenceGenerator_Pharr::minDist(SamplePMJ& pt, float* min_dist)
{
    float minSquareDist = 2; // On a unit grid. so can't be larger than this

    int c_i = int(grid_dim * pt[0] );
    int c_j = int(grid_dim * pt[1] );
    for (int i = -1; i <= 1; ++i)
    {
        int ii = c_i + i;
        if ( ii >= 0 && ii < grid_dim )
            for (int j = -1; j <= 1; ++j)
            {
                int jj = c_j + j;
                if ( jj >= 0 && jj < grid_dim )
                {
                    SamplePMJ* neighbor = grid_item( ii, jj );

                    if ( neighbor != nullptr )
                    {
                        float dist = pow((*neighbor)[0] - pt[0], 2) + 
                                      pow((*neighbor)[1] - pt[1], 2);
                        if(dist < minSquareDist)
                        {
                            minSquareDist = dist;
                        }
                    }
                }

            }
    }
    *min_dist = sqrt(minSquareDist);

    return ( minSquareDist > 1 ? false : true );
}

float PMJ02SampleSequenceGenerator_Pharr::generateRandomFloat() {
    return randomNumberGenerator->floatValue();
}

void PMJ02SampleSequenceGenerator_Pharr::initialize_x_tree( Node &node, int i, int j,
                                                            int shape, int nx,
                                                            int ny )
{

   node.min_val = i * ny;
   node.max_val = node.min_val + ny;

   if ( occupiedStrata[ shape ][ nx*j + i ] )
      node.occupied = true;
   else
   {
      node.occupied = false;
      node.leaf = ny == 1 ? true : false;
      if ( !node.leaf )
      {
         Node &node_left = x_tree[ ++x_tree_idx ];
         node.left = x_tree_idx;
         initialize_x_tree( node_left, 2 * i, int( j / 2 ), shape - 1, 2 * nx,
                            ny / 2 );

         Node &node_right = x_tree[ ++x_tree_idx ];
         node.right = x_tree_idx;
         initialize_x_tree( node_right, 2 * i + 1, int( j / 2 ), shape - 1, 2 * nx,
                            ny / 2 );
         node.occupied = node_left.occupied & node_right.occupied;
      }
   }
}

void PMJ02SampleSequenceGenerator_Pharr::initialize_y_tree( Node &node, int i, int j,
                                                            int shape, int nx,
                                                            int ny )
{

   node.min_val = j * nx;
   node.max_val = node.min_val + nx;
   if ( occupiedStrata[ shape ][ nx*j + i ] )
      node.occupied = true;
   else
   {
      node.occupied = false;
      node.leaf = nx == 1 ? true : false;
      if ( !node.leaf )
      {
         Node &node_left = y_tree[ ++y_tree_idx ];
         node.left = y_tree_idx;
         initialize_y_tree( node_left, int( i / 2 ), 2 * j, shape + 1, nx / 2,
                            2 * ny );

         Node &node_right = y_tree[ ++y_tree_idx ];
         node.right = y_tree_idx;
         initialize_y_tree( node_right, int( i / 2 ), 2 * j + 1, shape + 1, nx / 2,
                            2 * ny );
         node.occupied = node_left.occupied & node_right.occupied;
      }
   }
}

void
PMJ02SampleSequenceGenerator_Pharr::valid_x_stratum( Node& node, int i, int j, int nx, int ny )
{
    if ( !node.occupied )
    {
        if ( node.leaf )
            x_offsets.push_back(i);
        else
        {
            valid_x_stratum( x_tree[node.left],  2*i,     floor(j/2), 2*nx, ny/2 );
            valid_x_stratum( x_tree[node.right], 2*i + 1, floor(j/2), 2*nx, ny/2 );
        }
    }
}

void
PMJ02SampleSequenceGenerator_Pharr::valid_y_stratum( Node& node, int i, int j, int nx, int ny )
{
    if ( !node.occupied )
    {
        if ( node.leaf )
            y_offsets.push_back(j);
        else
        {
            valid_y_stratum( y_tree[node.left],  floor(i/2), 2*j,     nx/2, 2*ny );
            valid_y_stratum( y_tree[node.right], floor(i/2), 2*j + 1, nx/2, 2*ny );
        }
    }
}

//------------------------------------------------------------------------------
// PMJ02Sequence
//------------------------------------------------------------------------------

PMJ02Sequence::PMJ02Sequence(luxrays::RandomGenerator *rnd) : 
	rndGen(rnd), num_samples(4096) {
    localRng = false;

}

PMJ02Sequence::PMJ02Sequence(const u_int seed) : num_samples(4096) {
    localRng = true;
    rndGen = new luxrays::RandomGenerator(seed);
}

PMJ02Sequence::~PMJ02Sequence() {
    if (localRng) delete rndGen;
}

void PMJ02Sequence::RequestSamples(const u_int size) {
    boost::lock_guard<boost::mutex> guard(sampleGenerationMtx);
    // Samples already requested
    if (samplePoints.size() > 0) return; 
    
	// We cannot generate an odd number of dimensions
	u_int tablesToGenerate = (size / 2) + size % 2;

	samplePoints.resize(tablesToGenerate);

	SLG_LOG("Generating " << num_samples << " samples for " << size << " dimensions on " << tablesToGenerate << " tables");
	for (u_int i = 0; i < tablesToGenerate; i++) {
		samplePoints[i].resize(num_samples);
		PMJ02SampleSequenceGenerator_Pharr g = PMJ02SampleSequenceGenerator_Pharr(rndGen);
		g.ProgressiveMultiJittered02Algorithm2D(num_samples, 100);
		g.generatedSamples = shuffle(g.generatedSamples, num_samples);
		for (u_int j = 0; j < num_samples; j++) {
			SamplePMJ s = g.generatedSamples[j];
			samplePoints[i][j] = s;
		}
	}
	// SLG_LOG("after: " << samplePoints.size() << " " << samplePoints[0].size());
	for (u_int i = 0; i < tablesToGenerate; i++) {
		for (u_int j = 0; j < num_samples; j++) {
		}
	}
	SLG_LOG("Generated " << num_samples << " samples for " << size << " dimensions on " << tablesToGenerate << " tables");
}

float PMJ02Sequence::GetSample(const u_int pass, const u_int index) {

	// Revert to random if more samples than generated are being requested
	if (pass > num_samples) return rndGen->floatValue();

	const u_int dimensionIndex = (index / 2);
	return samplePoints[dimensionIndex][pass][index % 2];
}

std::vector<float> PMJ02Sequence::GetSamples(const u_int pass, const u_int offset) {

	std::vector<float> samples;
	samples.reserve(samplePoints.size() * 2);
	// TODO: Implement fallback after more samples than requested

    const u_int currentPass = (pass + offset) % samplePoints[0].size();
    // u_int currentPass = pass;
	for (u_int i = 0; i < samplePoints.size(); i++) {
		samples.push_back(samplePoints[i][currentPass][0]);
		samples.push_back(samplePoints[i][currentPass][1]);
	}

	return samples;
}

std::vector<SamplePMJ> PMJ02Sequence::shuffle(std::vector<SamplePMJ> points, u_int size) {

	constexpr u_int odd[8] = { 0, 1, 4, 5, 10, 11, 14, 15 };
	constexpr u_int even[8] = { 2, 3, 6, 7, 8, 9, 12, 13 };

	for (u_int yy = 0; yy < size / 16; ++yy) {
		for (u_int xx = 0; xx < 8; ++xx) {
			u_int other = (u_int)(rndGen->floatValue() * (8.0f - xx) + xx);
			SamplePMJ tmp = points[odd[other] + yy * 16];
			points[odd[other] + yy * 16] = points[odd[xx] + yy * 16];
			points[odd[xx] + yy * 16] = tmp;
		}
		for (u_int xx = 0; xx < 8; ++xx) {
			u_int other = (u_int)(rndGen->floatValue() * (8.0f - xx) + xx);
			SamplePMJ tmp = points[even[other] + yy * 16];
			points[even[other] + yy * 16] = points[even[xx] + yy * 16];
			points[even[xx] + yy * 16] = tmp;
		}
	}
    return points;
}

// std::vector<SamplePMJ> PMJ02Sequence::shuffle(std::vector<SamplePMJ> points, u_int size) {

//     for (u_int i = points.size(); i > 0; i--) { 
//         u_int j = rndGen->uintValue() % points.size();
//         const SamplePMJ tmp = points[i];
//         points[i] = points[j];
//         points[j] = tmp;
//     } 

//     return points;

// }
