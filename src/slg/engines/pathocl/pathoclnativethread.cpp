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

#if !defined(LUXRAYS_DISABLE_OPENCL)

#include <boost/lexical_cast.hpp>

#include "luxrays/core/geometry/transform.h"
#include "luxrays/core/randomgen.h"
#include "luxrays/utils/ocl.h"
#include "luxrays/core/oclintersectiondevice.h"
#include "luxrays/kernels/kernels.h"

#include "slg/slg.h"
#include "slg/kernels/kernels.h"
#include "slg/renderconfig.h"
#include "slg/engines/pathocl/pathocl.h"

using namespace std;
using namespace luxrays;
using namespace slg;

//------------------------------------------------------------------------------
// PathOCLNativeRenderThread
//------------------------------------------------------------------------------

PathOCLNativeRenderThread::PathOCLNativeRenderThread(const u_int index,
		NativeThreadIntersectionDevice *device, PathOCLRenderEngine *re) :
		PathOCLBaseNativeRenderThread(index, device, re) {
	threadFilm = NULL;
}

PathOCLNativeRenderThread::~PathOCLNativeRenderThread() {
	delete threadFilm; 
}

void PathOCLNativeRenderThread::Start() {
	if (threadIndex == 0) {
		// Only the first thread allocate a film. It is than used by all
		// other threads too.
		PathOCLRenderEngine *engine = (PathOCLRenderEngine *)renderEngine;

		const u_int filmWidth = engine->film->GetWidth();
		const u_int filmHeight = engine->film->GetHeight();
		const u_int *filmSubRegion = engine->film->GetSubRegion();

		delete threadFilm;

		threadFilm = new Film(filmWidth, filmHeight, filmSubRegion);
		threadFilm->CopyDynamicSettings(*(engine->film));
		// I'm not removing the pipeline and disabling the film denoiser
		// in order to support BCD denoiser.
		threadFilm->Init();
	}

	PathOCLBaseNativeRenderThread::Start();
}

void PathOCLNativeRenderThread::StartRenderThread() {
	if (threadIndex == 0) {
		// Clear the film (i.e. the thread can be started/stopped multiple times)
		threadFilm->Clear();
	}

	PathOCLBaseNativeRenderThread::StartRenderThread();
}

void PathOCLNativeRenderThread::RenderThreadImpl() {
	//SLG_LOG("[PathOCLRenderEngine::" << threadIndex << "] Rendering thread started");

	// Boost barriers (used in PhotonGICache::Update()) are supposed to be not
	// interruptible but they are and seem to be missing a way to reset them. So
	// better to disable interruptions.
	boost::this_thread::disable_interruption di;

	//--------------------------------------------------------------------------
	// Initialization
	//--------------------------------------------------------------------------

	PathOCLRenderEngine *engine = (PathOCLRenderEngine *)renderEngine;
	const PathTracer &pathTracer = engine->pathTracer;
	// (engine->seedBase + 1) seed is used for sharedRndGen
	RandomGenerator *rndGen = new RandomGenerator(engine->seedBase + 1 + threadIndex);

	// All threads use the film allocated by the first thread
	Film *film = ((PathOCLNativeRenderThread *)(engine->renderNativeThreads[0]))->threadFilm;
	
	// Setup the sampler(s)

	Sampler *eyeSampler = nullptr;
	Sampler *lightSampler = nullptr;

	eyeSampler = engine->renderConfig->AllocSampler(rndGen, film,
			nullptr, engine->eyeSamplerSharedData, Properties());
	eyeSampler->RequestSamples(PIXEL_NORMALIZED_ONLY, pathTracer.eyeSampleSize);

	if (pathTracer.hybridBackForwardEnable) {
		// Light path sampler is always Metropolis
		Properties props;
		props <<
			Property("sampler.type")("METROPOLIS") <<
			// Disable image plane meaning for samples 0 and 1
			Property("sampler.imagesamples.enable")(false);

		lightSampler = Sampler::FromProperties(props, rndGen, film, nullptr,
				engine->lightSamplerSharedData);
		
		lightSampler->RequestSamples(SCREEN_NORMALIZED_ONLY, pathTracer.lightSampleSize);
	}
	
	VarianceClamping varianceClamping(pathTracer.sqrtVarianceClampMaxValue);

	// Setup PathTracer thread state
	PathTracerThreadState pathTracerThreadState(threadIndex, intersectionDevice,
			eyeSampler, lightSampler,
			engine->renderConfig->scene, film,
			&varianceClamping);

	//--------------------------------------------------------------------------
	// Trace paths
	//--------------------------------------------------------------------------

	// I can not use engine->renderConfig->GetProperty() here because the
	// RenderConfig properties cache is not thread safe
	const u_int filmWidth = film->GetWidth();
	const u_int filmHeight = film->GetHeight();
	const u_int haltDebug = engine->renderConfig->cfg.Get(Property("batch.haltdebug")(0u)).Get<u_int>() *
		filmWidth * filmHeight;

	for (u_int steps = 0; !boost::this_thread::interruption_requested(); ++steps) {
		// Check if we are in pause mode
		if (engine->pauseMode) {
			// Check every 100ms if I have to continue the rendering
			while (!boost::this_thread::interruption_requested() && engine->pauseMode)
				boost::this_thread::sleep(boost::posix_time::millisec(100));

			if (boost::this_thread::interruption_requested())
				break;
		}

		pathTracer.RenderSample(pathTracerThreadState);

#ifdef WIN32
		// Work around Windows bad scheduling
		renderThread->yield();
#endif

		// Check halt conditions
		if ((haltDebug > 0u) && (steps >= haltDebug))
			break;
		if (engine->film->GetConvergence() == 1.f)
			break;

		if (engine->photonGICache)
			engine->photonGICache->Update(engine->renderOCLThreads.size() + threadIndex, engine->GetTotalEyeSPP());
	}

	delete eyeSampler;
	delete lightSampler;
	delete rndGen;

	threadDone = true;

	//SLG_LOG("[PathOCLRenderEngine::" << threadIndex << "] Rendering thread halted");
}

#endif
