#line 2 "materialdefs_funcs_glass.cl"

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

//------------------------------------------------------------------------------
// Glass material
//------------------------------------------------------------------------------

#if defined (PARAM_ENABLE_MAT_GLASS)

OPENCL_FORCE_INLINE BSDFEvent GlassMaterial_GetEventTypes() {
	return SPECULAR | REFLECT | TRANSMIT;
}

OPENCL_FORCE_INLINE bool GlassMaterial_IsDelta() {
	return true;
}

OPENCL_FORCE_INLINE float3 GlassMaterial_Evaluate(
		__global const HitPoint *hitPoint, const float3 lightDir, const float3 eyeDir,
		BSDFEvent *event, float *directPdfW,
		const float3 ktTexVal, const float3 krTexVal,
		const float3 nc, const float3 nt, const float cauchyC) {
	return BLACK;
}

OPENCL_FORCE_INLINE float3 GlassMaterial_WaveLength2RGB(const float waveLength) {
	float r, g, b;
	if ((waveLength >= 380.f) && (waveLength < 440.f)) {
		r = -(waveLength - 440.f) / (440 - 380.f);
		g = 0.f;
		b = 1.f;
	} else if ((waveLength >= 440.f) && (waveLength < 490.f)) {
		r = 0.f;
		g = (waveLength - 440.f) / (490.f - 440.f);
		b = 1.f;
	} else if ((waveLength >= 490.f) && (waveLength < 510.f)) {
		r = 0.f;
		g = 1.f;
		b = -(waveLength - 510.f) / (510.f - 490.f);
	} else if ((waveLength >= 510.f) && (waveLength < 580.f)) {
		r = (waveLength - 510.f) / (580.f - 510.f);
		g = 1.f;
		b = 0.f;
	} else if ((waveLength >= 580.f) && (waveLength < 645.f)) {
		r = 1.f;
		g = -(waveLength - 645.f) / (645 - 580.f);
		b = 0.f;
	} else if ((waveLength >= 645.f) && (waveLength < 780.f)) {
		r = 1.f;
		g = 0.f;
		b = 0.f;
	} else
		return BLACK;

	// The intensity fall off near the upper and lower limits
	float factor;
	if ((waveLength >= 380.f) && (waveLength < 420.f))
		factor = .3f + .7f * (waveLength - 380.f) / (420.f - 380.f);
	else if ((waveLength >= 420) && (waveLength < 700))
		factor = 1.f;
	else
		factor = .3f + .7f * (780.f - waveLength) / (780.f - 700.f);

	const float3 result = (float3)(r, g, b) * factor;

	/*
	Spectrum white;
	for (u_int i = 380; i < 780; ++i)
		white += WaveLength2RGB(i);
	white *= 1.f / 400.f;
	cout << std::setprecision(std::numeric_limits<float>::digits10 + 1) << white.c[0] << ", " << white.c[1] << ", " << white.c[2] << "\n";
	 
	 Result: 0.5652729, 0.36875, 0.265375
	 */

	// To normalize the output
	const float3 normFactor = (float3)(1.f / .5652729f, 1.f / .36875f, 1.f / .265375f);
	
	return result * normFactor;
}

OPENCL_FORCE_INLINE float GlassMaterial_WaveLength2IOR(const float waveLength, const float IOR, const float C) {
	// Cauchy's equation for relationship between the refractive index and wavelength
	// note: Cauchy's lambda is expressed in micrometers while waveLength is in nanometers

	// This is the formula suggested by Neo here:
	// https://github.com/LuxCoreRender/BlendLuxCore/commit/d3fed046ab62e18226e410b42a16ca1bccefb530#commitcomment-26617643
	//const float B = IOR - C / Sqr(589.f / 1000.f);

	// The B used by old LuxRender
	const float B = IOR;

	// Cauchy's equation
	const float cauchyEq = B + C / Sqr(waveLength / 1000.f);

	return cauchyEq;
}

OPENCL_FORCE_INLINE float3 GlassMaterial_EvalSpecularReflection(__global const HitPoint *hitPoint,
		const float3 localFixedDir, const float3 kr,
		const float nc, const float nt,
		float3 *localSampledDir) {
	if (Spectrum_IsBlack(kr))
		return BLACK;

	const float costheta = CosTheta(localFixedDir);
	*localSampledDir = (float3)(-localFixedDir.x, -localFixedDir.y, localFixedDir.z);

	const float ntc = nt / nc;
	return kr * FresnelCauchy_Evaluate(ntc, costheta);
}

OPENCL_FORCE_INLINE float3 GlassMaterial_EvalSpecularTransmission(__global const HitPoint *hitPoint,
		const float3 localFixedDir, const float u0,
		const float3 kt, const float nc, const float nt, const float cauchyC,
		float3 *localSampledDir) {
	if (Spectrum_IsBlack(kt))
		return BLACK;

	// Compute transmitted ray direction
	float3 lkt;
	float lnt;
	if (cauchyC > 0.f) {
		// Select the wavelength to sample
		const float waveLength = mix(380.f, 780.f, u0);

		lnt = GlassMaterial_WaveLength2IOR(waveLength, nt, cauchyC);

		lkt = kt * GlassMaterial_WaveLength2RGB(waveLength);
	} else {
		lnt = nt;
		lkt = kt;
	}

	const float ntc = lnt / nc;
	const float costheta = CosTheta(localFixedDir);
	const bool entering = (costheta > 0.f);
	const float eta = entering ? (nc / lnt) : ntc;
	const float eta2 = eta * eta;
	const float sini2 = SinTheta2(localFixedDir);
	const float sint2 = eta2 * sini2;

	// Handle total internal reflection for transmission
	if (sint2 >= 1.f)
		return BLACK;

	const float cost = sqrt(fmax(0.f, 1.f - sint2)) * (entering ? -1.f : 1.f);
	*localSampledDir = (float3)(-eta * localFixedDir.x, -eta * localFixedDir.y, cost);

	float ce;
//	if (!hitPoint.fromLight)
		ce = (1.f - FresnelCauchy_Evaluate(ntc, cost)) * eta2;
//	else {
//		const float absCosSampledDir = fabsf(CosTheta(*localSampledDir));
//		ce = (1.f - FresnelTexture::CauchyEvaluate(ntc, costheta)) * fabsf(localFixedDir.z / absCosSampledDir);
//	}

	return lkt * ce;
}

OPENCL_FORCE_NOT_INLINE float3 GlassMaterial_Sample(
		__global const HitPoint *hitPoint, const float3 localFixedDir, float3 *localSampledDir,
		const float u0, const float u1,
		const float passThroughEvent,
		float *pdfW, BSDFEvent *event,
		const float3 ktTexVal, const float3 krTexVal,
		const float nc, const float nt, const float cauchyC) {
	const float3 kt = Spectrum_Clamp(ktTexVal);
	const float3 kr = Spectrum_Clamp(krTexVal);

	float3 transLocalSampledDir; 
	const float3 trans = GlassMaterial_EvalSpecularTransmission(hitPoint, localFixedDir, u0,
			kt, nc, nt, cauchyC, &transLocalSampledDir);
	
	float3 reflLocalSampledDir;
	const float3 refl = GlassMaterial_EvalSpecularReflection(hitPoint, localFixedDir,
			kr, nc, nt, &reflLocalSampledDir);

	// Decide to transmit or reflect
	float threshold;
	if (!Spectrum_IsBlack(refl)) {
		if (!Spectrum_IsBlack(trans)) {
			// Importance sampling
			const float reflFilter = Spectrum_Filter(refl);
			const float transFilter = Spectrum_Filter(trans);
			threshold = transFilter / (reflFilter + transFilter);
		} else
			threshold = 0.f;
	} else {
		if (!Spectrum_IsBlack(trans))
			threshold = 1.f;
		else
			return BLACK;
	}

	float3 result;
	if (passThroughEvent < threshold) {
		// Transmit

		*localSampledDir = transLocalSampledDir;

		*event = SPECULAR | TRANSMIT;
		*pdfW = threshold;
	
		result = trans;
	} else {
		// Reflect

		*localSampledDir = reflLocalSampledDir;

		*event = SPECULAR | REFLECT;
		*pdfW = 1.f - threshold;
		
		result = refl;
	}

	return result / *pdfW;
}

#endif
