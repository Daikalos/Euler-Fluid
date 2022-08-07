#include "Camera.h"

Camera::Camera(const sf::RenderWindow& window) : window(window), position(sf::Vector2f(window.getSize()) / 2.0f), scale({ 1.0f, 1.0f })
{

}

void Camera::update(const InputHandler& inputHandler)
{
	if (inputHandler.get_key_pressed(sf::Keyboard::Key::Space))
	{
		position = (sf::Vector2f)window.getSize() / 2.0f;
	}

	if (inputHandler.get_middle_pressed())
		dragPos = get_mouse_world_position();
	if (inputHandler.get_middle_held())
		position += (sf::Vector2f)(dragPos - get_mouse_world_position());
}
