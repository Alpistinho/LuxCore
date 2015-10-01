/***************************************************************************
 * Copyright 1998-2015 by authors (see AUTHORS.txt)                        *
 *                                                                         *
 *   This file is part of LuxRender.                                       *
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

#ifndef _SLG_SOBOL_SAMPLER_H
#define	_SLG_SOBOL_SAMPLER_H

#include <string>
#include <vector>

#include "luxrays/core/randomgen.h"
#include "slg/slg.h"
#include "slg/film/film.h"
#include "slg/samplers/sampler.h"

namespace slg {

//------------------------------------------------------------------------------
// SobolSamplerSharedData
//
// Used to share sampler specific data across multiple threads
//------------------------------------------------------------------------------

class SobolSamplerSharedData : public SamplerSharedData{
public:
	SobolSamplerSharedData(luxrays::RandomGenerator *rnd);
	virtual ~SobolSamplerSharedData() { }

	float rng0, rng1;
};

//------------------------------------------------------------------------------
// Sobol sampler
//
// This sampler is based on Blender Cycles Sobol implementation.
//------------------------------------------------------------------------------

#define SOBOL_STARTOFFSET 32

extern void SobolGenerateDirectionVectors(u_int *vectors, const u_int dimensions);

class SobolSampler : public Sampler {
public:
	SobolSampler(luxrays::RandomGenerator *rnd, Film *flm,
			const u_int threadIndex, const u_int threadCount,
			SobolSamplerSharedData *samplerSharedData);
	virtual ~SobolSampler();

	virtual SamplerType GetType() const { return SOBOL; }
	virtual void RequestSamples(const u_int size);

	virtual float GetSample(const u_int index);
	virtual void NextSample(const std::vector<SampleResult> &sampleResults);

	static SobolSamplerSharedData *AllocSharedData(luxrays::RandomGenerator *rnd);

private:
	u_int SobolDimension(const u_int index, const u_int dimension) const;

	SobolSamplerSharedData *sharedData;
	const u_int threadIndex, threadCount;

	u_int *directions;
	u_int pass;
};

}

#endif	/* _SLG_SOBOL_SAMPLER_H */
