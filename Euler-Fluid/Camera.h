#pragma once

#include <SFML/Graphics.hpp>
#include <SFML/OpenGL.hpp>
#include <gl/GLU.h>
#include <iostream>

#include "InputHandler.h"

class Camera
{
public:
	Camera() = delete;
	Camera(const sf::Window& window);

	// call after poll event
	//
	void update(const InputHandler& inputHandler);

	template<typename T> sf::Vector2<T> view_to_world(const sf::Vector2<T>& position) const
	{
		return (sf::Vector2<T>)(get_view_matrix() * (sf::Vector2f)position);
	}

public:
	inline sf::View get_view() const
	{
		return view;
	}

	inline const float* get_world_matrix() const
	{
		return sf::Transform()
			.translate((sf::Vector2f)window.getSize() / 2.0f)
			.scale(scale, scale)
			.translate(-position).getMatrix();
	}
	inline sf::Transform get_view_matrix() const
	{
		return sf::Transform()
			.translate(position)
			.scale(1.0f / scale, 1.0f / scale)
			.translate((sf::Vector2f)window.getSize() / -2.0f);
	}

	inline sf::Vector2f get_position() const { return position; }
	inline sf::Vector2i get_mouse_world_position() const { return view_to_world(sf::Mouse::getPosition(window)); }

	inline double get_scale() const { return scale; }

private:
	const sf::Window& window;

	sf::Vector2f position;
	sf::Vector2f size;
	float scale;

	sf::Vector2i dragPos;

	sf::View view;
};

