#pragma once

#define _USE_MATH_DEFINES

#include <math.h>

namespace util
{
	template<typename T> static inline T to_radians(T degrees)
	{
		return (T)(degrees * ((double)M_PI / 180.0));
	}

	template<typename T> static inline T to_degrees(T radians)
	{
		return (T)(radians * (180.0 / (double)M_PI));
	}

	template<typename T> static inline T clamp(T val, const T min, const T max)
	{
		val = (val < min) ? min : val;
		val = (val > max) ? max : val;
		return val;
	}

	template<typename T> static inline T fRand(T fMin, T fMax)
	{
		T f = (T)rand() / RAND_MAX;
		return fMin + f * (fMax - fMin);
	}
}
