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
#include <string>

#include <boost/lexical_cast.hpp>

#include "luxrays/core/color/color.h"
#include "slg/samplers/sampler.h"
#include "slg/samplers/pmj02.h"

using namespace std;
using namespace luxrays;
using namespace slg;

u_int cmj_hash_simple(u_int i, u_int p) {
	i = (i ^ 61) ^ p;
	i += i << 3;
	i ^= i >> 4;
	i *= 0x27d4eb2d;
	return i;
}

inline uint64_t rotl(const uint64_t x, int k) {
	return (x << k) | (x >> (64 - k));
}


//------------------------------------------------------------------------------
// PMJ02SamplerSharedData
//------------------------------------------------------------------------------


PMJ02SamplerSharedData::PMJ02SamplerSharedData(const u_int seed, Film *engineFlm) : SamplerSharedData(), pmj02sequence(seed) {
	Init(seed, engineFlm);
}

PMJ02SamplerSharedData::PMJ02SamplerSharedData(RandomGenerator *rndGen, Film *engineFlm) : SamplerSharedData(), pmj02sequence(rndGen) {
	Init(rndGen->uintValue() % (0xFFFFFFFFu - 1u) + 1u, engineFlm);
}

void PMJ02SamplerSharedData::Init(const u_int seed, Film *engineFlm) {
	engineFilm = engineFlm;
	seedBase = seed;

	if (engineFilm) {
		const u_int *subRegion = engineFilm->GetSubRegion();
		filmRegionPixelCount = (subRegion[1] - subRegion[0] + 1) * (subRegion[3] - subRegion[2] + 1);

		// All passes start at 0
		passPerPixel.resize(filmRegionPixelCount, 0);
	} else {
		filmRegionPixelCount = 0;
		passPerPixel.resize(1, 0);
	}

	pixelIndex = 0;
}

void PMJ02SamplerSharedData::GetNewPixelIndex(u_int &index, u_int &seed) {
	SpinLocker spinLocker(spinLock);

	index = pixelIndex;
	seed = (seedBase + pixelIndex) % (0xFFFFFFFFu - 1u) + 1u;

	pixelIndex += PMJ02_THREAD_WORK_SIZE;
	if (pixelIndex >= filmRegionPixelCount)
		pixelIndex = 0;
}

u_int PMJ02SamplerSharedData::GetNewPixelPass(const u_int pixelIndex) {
	// Iterate pass of this pixel
	return AtomicInc(&passPerPixel[pixelIndex]);
}

float PMJ02SamplerSharedData::GetSample(const u_int pass, const u_int index) {
	// boost::shared_lock<boost::shared_mutex> samplesReady(sampleGenerationMutex);
	return pmj02sequence.GetSample(pass, index);;
}

std::vector<float> PMJ02SamplerSharedData::GetSamples(const u_int pass, const u_int offset) {
	return pmj02sequence.GetSamples(pass, offset);
}

void PMJ02SamplerSharedData::RequestSamples(const u_int size) {
	boost::unique_lock<boost::shared_mutex> samplesReady(sampleGenerationMutex);

	requestedSamples = size;
    pmj02sequence.RequestSamples(size);
}

SamplerSharedData *PMJ02SamplerSharedData::FromProperties(const Properties &cfg,
		RandomGenerator *rndGen, Film *film) {
	return new PMJ02SamplerSharedData(rndGen, film);
}



//------------------------------------------------------------------------------
// PMJ02 sampler
//------------------------------------------------------------------------------

PMJ02Sampler::PMJ02Sampler(luxrays::RandomGenerator *rnd, Film *flm,
			const FilmSampleSplatter *flmSplatter, const bool imgSamplesEnable,
			const float adaptiveStr,
			PMJ02SamplerSharedData *samplerSharedData) :
		Sampler(rnd, flm, flmSplatter, imgSamplesEnable), sharedData(samplerSharedData),
		adaptiveStrength(adaptiveStr) {
}

void PMJ02Sampler::InitNewSample() {
	for (;;) {
		// Update pixelIndexOffset

		pixelIndexOffset++;
		if (pixelIndexOffset >= PMJ02_THREAD_WORK_SIZE) {
			// Ask for a new base
			u_int seed;
			sharedData->GetNewPixelIndex(pixelIndexBase, seed);
			pixelIndexOffset = 0;
		}

		// Initialize sample0 and sample 1

		u_int pixelX, pixelY;
		if (imageSamplesEnable && film) {
			const u_int *subRegion = film->GetSubRegion();

			pixelIndex = (pixelIndexBase + pixelIndexOffset) % sharedData->filmRegionPixelCount;
			const u_int subRegionWidth = subRegion[1] - subRegion[0] + 1;
			pixelX = subRegion[0] + (pixelIndex % subRegionWidth);
			pixelY = subRegion[2] + (pixelIndex / subRegionWidth);

			// Check if the current pixel is over or under the noise threshold
			const Film *film = sharedData->engineFilm;
			if ((adaptiveStrength > 0.f) && film->HasChannel(Film::NOISE)) {
				// Pixels are sampled in accordance with how far from convergence they are
				// The floor for the pixel importance is given by the adaptiveness strength
				const float noise = Max(*(film->channel_NOISE->GetPixel(pixelX, pixelY)), 1.f - adaptiveStrength);
				pass = sharedData->GetNewPixelPass(pixelIndex);
				if (rndGen->floatValue() > noise) {
					// Skip this pixel and try the next one
					pass = sharedData->GetNewPixelPass(pixelIndex);
					continue;
				}
			}

		} else {
			pixelX = 0;
			pixelY = 0;
		}
		seed[0] = pixelIndex;
		seed[1] = 534613516;
		
		currentSamples = sharedData->GetSamples(pass, next_random());
		
		// Shuffle array using Peter-Yates algorithm
		for (u_int i = (currentSamples.size()/2) - 1; i > 0; i--) { 
			u_int j = next_random() % (currentSamples.size()/2);
			j = j - (j % 2);

			const u_int currentDimension = 2*i;
			const u_int swapDimension = 2*j;
			const float tmp1 = currentSamples[currentDimension];
			const float tmp2 = currentSamples[currentDimension + 1];
			currentSamples[currentDimension] = currentSamples[swapDimension];
			currentSamples[currentDimension + 1] = currentSamples[swapDimension + 1];
			currentSamples[swapDimension] = tmp1;
			currentSamples[swapDimension + 1] = tmp2;
		} 

		// if (pixelIndex == 12312) {
		// 	std::ostringstream ss;
		// 	for (u_int i = 0; i < currentSamples.size(); i++ ) {
		// 		ss << currentSamples[i] << ",";
		// 	}
		// 	std::string s(ss.str());
		// 	SLG_LOG("Sample: " << s);
		// 	SLG_LOG("Next_random: " << next_random());
		// } 

		sample0 = pixelX + currentSamples[0];
		sample1 = pixelY + currentSamples[1];
		break;
	}
}

void PMJ02Sampler::RequestSamples(const SampleType smplType, const u_int size) {
	Sampler::RequestSamples(smplType, size);
	sharedData->RequestSamples(size);
	pixelIndexOffset = PMJ02_THREAD_WORK_SIZE;
	InitNewSample();
}

float PMJ02Sampler::GetSample(const u_int index) {
	assert (index < requestedSamples);
	switch (index) {
		case 0:
			return sample0;
		case 1:
			return sample1;
		default:
			return currentSamples[index];
	}
}

void PMJ02Sampler::NextSample(const vector<SampleResult> &sampleResults) {
	if (film) {
		double pixelNormalizedCount, screenNormalizedCount;
		switch (sampleType) {
			case PIXEL_NORMALIZED_ONLY:
				pixelNormalizedCount = 1.0;
				screenNormalizedCount = 0.0;
				break;
			case SCREEN_NORMALIZED_ONLY:
				pixelNormalizedCount = 0.0;
				screenNormalizedCount = 1.0;
				break;
			case PIXEL_NORMALIZED_AND_SCREEN_NORMALIZED:
				pixelNormalizedCount = 1.0;
				screenNormalizedCount = 1.0;
				break;
			default:
				throw runtime_error("Unknown sample type in PMJ02Sampler::NextSample(): " + ToString(sampleType));
		}
		film->AddSampleCount(threadIndex, pixelNormalizedCount, screenNormalizedCount);

		AtomicAddSamplesToFilm(sampleResults);
	}

	InitNewSample();
}

inline uint64_t PMJ02Sampler::next_random(void) {
	const uint64_t s0 = seed[0];
	uint64_t s1 = seed[1];
	const uint64_t result = rotl(s0 * 5, 7) * 9;

	s1 ^= s0;
	seed[0] = rotl(s0, 24) ^ s1 ^ (s1 << 16); // a, b
	seed[1] = rotl(s1, 37); // c

	return result;
}

Properties PMJ02Sampler::ToProperties() const {
	return Sampler::ToProperties() <<
			Property("sampler.pmj02.adaptive.strength")(adaptiveStrength);
}

//------------------------------------------------------------------------------
// Static methods used by SamplerRegistry
//------------------------------------------------------------------------------

Properties PMJ02Sampler::ToProperties(const Properties &cfg) {
	return Properties() <<
			cfg.Get(GetDefaultProps().Get("sampler.type")) <<
			cfg.Get(GetDefaultProps().Get("sampler.imagesamples.enable")) <<
			cfg.Get(GetDefaultProps().Get("sampler.pmj02.adaptive.strength"));
}

Sampler *PMJ02Sampler::FromProperties(const Properties &cfg, RandomGenerator *rndGen,
		Film *film, const FilmSampleSplatter *flmSplatter, SamplerSharedData *sharedData) {
	const bool imageSamplesEnable = cfg.Get(GetDefaultProps().Get("sampler.imagesamples.enable")).Get<bool>();

	const float str = Clamp(cfg.Get(GetDefaultProps().Get("sampler.pmj02.adaptive.strength")).Get<float>(), 0.f, .95f);

	return new PMJ02Sampler(rndGen, film, flmSplatter, imageSamplesEnable,
			str, (PMJ02SamplerSharedData *)sharedData);
}

slg::ocl::Sampler *PMJ02Sampler::FromPropertiesOCL(const Properties &cfg) {
	slg::ocl::Sampler *oclSampler = new slg::ocl::Sampler();

	oclSampler->type = slg::ocl::PMJ02;
	oclSampler->pmj02.adaptiveStrength = Clamp(cfg.Get(GetDefaultProps().Get("sampler.pmj02.adaptive.strength")).Get<float>(), 0.f, .95f);

	return oclSampler;
}

Film::FilmChannelType PMJ02Sampler::GetRequiredChannels(const luxrays::Properties &cfg) {
	const bool imageSamplesEnable = cfg.Get(GetDefaultProps().Get("sampler.imagesamples.enable")).Get<bool>();

	const float str = cfg.Get(GetDefaultProps().Get("sampler.pmj02.adaptive.strength")).Get<float>();

	if (imageSamplesEnable && (str > 0.f))
		return Film::NOISE;
	else
		return Film::NONE;
}

const Properties &PMJ02Sampler::GetDefaultProps() {
	static Properties props = Properties() <<
			Sampler::GetDefaultProps() <<
			Property("sampler.type")(GetObjectTag()) <<
			Property("sampler.pmj02.adaptive.strength")(.95f);

	return props;
}