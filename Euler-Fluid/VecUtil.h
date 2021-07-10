#pragma once

#define _USE_MATH_DEFINES

#include <SFML/Graphics.hpp>
#include "Utilities.h"

template<typename T> struct v2
{
public:
	static inline float length(const sf::Vector2<T>& vector)
	{
		return std::sqrt(vector.x * vector.x + vector.y * vector.y);
	}
	static inline float angle(const sf::Vector2<T>& vector)
	{
		return atan2f(vector.y, vector.x);
	}
	static inline float angle(const sf::Vector2<T>& lhs, const sf::Vector2<T>& rhs)
	{
		return acosf(util::clamp(dot(lhs, rhs) / (length(lhs) * length(rhs)), -1.0f, 1.0f));
	}
	static inline float dot(const sf::Vector2<T>& lhs, const sf::Vector2<T>& rhs)
	{
		return lhs.x * rhs.x + lhs.y * rhs.y;
	}
	static inline float distance(const sf::Vector2<T>& from, const sf::Vector2<T>& to)
	{
		return length(direction(from, to));
	}

	static inline sf::Vector2<T> normalize(sf::Vector2<T> vector, float radius = 1.0f)
	{
		float len = length(vector);

		if (len < FLT_EPSILON)
			return vector;

		vector.x *= (radius / len);
		vector.y *= (radius / len);

		return vector;
	}
	static inline sf::Vector2<T> limit(const sf::Vector2<T>& vector, float maxLength)
	{
		if (length(vector) > maxLength)
		{
			return normalize(vector, maxLength);
		}
		return vector;
	}
	static inline sf::Vector2<T> direction(const sf::Vector2<T>& from, const sf::Vector2<T>& to)
	{
		return sf::Vector2<T>(to.x - from.x, to.y - from.y);
	}
	static inline sf::Vector2<T> rotate_point(const sf::Vector2<T>& point, const sf::Vector2<T>& center, float angle)
	{
		sf::Vector2<T> newPoint = sf::Vector2<T>();

		float s = sinf(angle);
		float c = cosf(angle);

		newPoint.x = (((point.x - center.x) * c - (point.y - center.y) * s)) + center.x;
		newPoint.y = (((point.x - center.x) * s + (point.y - center.y) * c)) + center.y;

		return newPoint;
	}

private:
	v2() = delete;
};

typedef v2<int> v2i;
typedef v2<float> v2f;
typedef v2<double> v2d;

template <typename T> 
static sf::Vector2<T> operator /=(sf::Vector2<T>& lhs, const sf::Vector2<T>& rhs)
{
	lhs.x /= rhs.x;
	lhs.y /= rhs.y;
	return lhs;
}

template <typename T>
static sf::Vector2<T> operator /(const sf::Vector2<T>& lhs, const sf::Vector2<T>& rhs)
{
	sf::Vector2<T> result = lhs;
	result.x /= rhs.x;
	result.y /= rhs.y;
	return result;
}

