#pragma once

#define _USE_MATH_DEFINES

#include <SFML/Graphics.hpp>

#include "Utilities.h"

template<typename T> struct v2
{
public:
	static constexpr float length(const sf::Vector2<T>& vector)
	{
		return std::sqrtf(vector.x * vector.x + vector.y * vector.y);
	}
	static constexpr float angle(const sf::Vector2<T>& vector)
	{
		return atan2f(vector.y, vector.x);
	}
	static constexpr float angle(const sf::Vector2<T>& lhs, const sf::Vector2<T>& rhs)
	{
		return acosf(util::clamp(dot(lhs, rhs) / (length(lhs) * length(rhs)), -1.0f, 1.0f));
	}
	static constexpr float dot(const sf::Vector2<T>& lhs, const sf::Vector2<T>& rhs)
	{
		return lhs.x * rhs.x + lhs.y * rhs.y;
	}
	static constexpr float distance(const sf::Vector2<T>& from, const sf::Vector2<T>& to)
	{
		T dx = std::abs(to.x - from.x);
		T dy = std::abs(to.y - from.y);

		if (dy > dx)
			std::swap(dx, dy);

		return 1007.f / 1024.f * dx + 441.f / 1024.f * dy;
	}
	static constexpr float distance_squared(const sf::Vector2<T>& from, const sf::Vector2<T>& to)
	{
		sf::Vector2f dir = direction(from, to);
		return dir.x * dir.x + dir.y * dir.y;
	}

	static constexpr sf::Vector2<T> normalize(sf::Vector2<T> vector, const float& radius = 1.0f)
	{
		float len = length(vector);

		if (len < FLT_EPSILON)
			return vector;

		vector.x *= (radius / len);
		vector.y *= (radius / len);

		return vector;
	}
	static constexpr sf::Vector2<T> limit(const sf::Vector2<T>& vector, const float& maxLength)
	{
		if (length(vector) > maxLength)
			return normalize(vector, maxLength);

		return vector;
	}
	static constexpr sf::Vector2<T> min(const sf::Vector2<T>& vector, const float& minLength)
	{
		if (length(vector) < minLength)
			return normalize(vector, minLength);

		return vector;
	}
	static constexpr sf::Vector2<T> direction(const sf::Vector2<T>& from, const sf::Vector2<T>& to)
	{
		return sf::Vector2<T>(to.x - from.x, to.y - from.y);
	}
	static constexpr sf::Vector2<T> rotate_point(const sf::Vector2<T>& point, const sf::Vector2<T>& center, const float& angle)
	{
		sf::Vector2<T> dir = direction(center, point);

		float s = sinf(angle);
		float c = cosf(angle);

		return sf::Vector2<T>(
			(dir.x * c - dir.y * s) + center.x,
			(dir.x * s + dir.y * c) + center.y);
	}
	static constexpr float lerp(float a, float b, float f)
	{
		return (a * (1.0 - f)) + (b * f);
	}
	static constexpr sf::Vector2<T> lerp(sf::Vector2<T> lhs, sf::Vector2<T> rhs, float a)
	{
		return sf::Vector3<T>(
			lerp(lhs.x, rhs.x, a),
			lerp(lhs.y, rhs.y, a));
	}
	static constexpr sf::Vector3<T> lerp(sf::Vector3<T> lhs, sf::Vector3<T> rhs, float a)
	{
		return sf::Vector3<T>(
			lerp(lhs.x, rhs.x, a),
			lerp(lhs.y, rhs.y, a),
			lerp(lhs.z, rhs.z, a));
	}

private:
	v2() = delete;
};

typedef v2<int> v2i;
typedef v2<float> v2f;
typedef v2<double> v2d;

template <typename T> 
static constexpr sf::Vector2<T> operator /=(sf::Vector2<T>& lhs, const sf::Vector2<T>& rhs)
{
	lhs.x /= rhs.x;
	lhs.y /= rhs.y;
	return lhs;
}
template <typename T>
static constexpr sf::Vector2<T> operator /(float lhs, const sf::Vector2<T>& rhs)
{
	sf::Vector2<T> result = { lhs, lhs };
	result.x /= rhs.x;
	result.y /= rhs.y;
	return result;
}

template <typename T>
static constexpr sf::Vector2<T> operator /(sf::Vector2<T> lhs, const sf::Vector2<T>& rhs)
{
	lhs.x /= rhs.x;
	lhs.y /= rhs.y;
	return lhs;
}

template <typename T>
static constexpr sf::Vector2<T> operator *(sf::Vector2<T> lhs, const sf::Vector2<T>& rhs)
{
	lhs.x *= rhs.x;
	lhs.y *= rhs.y;
	return lhs;
}
