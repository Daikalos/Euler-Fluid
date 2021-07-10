#include "Camera.h"

Camera::Camera(const sf::Window& window) : window(window)
{
	position = (sf::Vector2f)window.getSize() / 2.0f;
	scale = 1.0f;

	view.setCenter(position);
	view.setSize(size);
	view.setViewport(sf::FloatRect(0.0f, 0.0f, 1.0f, 1.0f));
}

void Camera::update(const InputHandler& inputHandler)
{
	if (inputHandler.key_pressed(sf::Keyboard::Key::Space))
	{
		position = (sf::Vector2f)window.getSize() / 2.0f;
		scale = 1.0f;
	}

	if (inputHandler.middle_pressed())
		dragPos = get_mouse_world_position();
	if (inputHandler.middle_held())
		position += (sf::Vector2f)(dragPos - get_mouse_world_position());

	if (inputHandler.scroll_down())
		scale *= 0.85f;
	if (inputHandler.scroll_up())
		scale *= 1.15f;

	view.setCenter(position);
	view.setSize((sf::Vector2f)window.getSize() * (1.0f / scale));
}
