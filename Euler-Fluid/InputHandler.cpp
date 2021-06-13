#include "InputHandler.h"

InputHandler::InputHandler()
{
	memset(currentButtonState, false, sizeof(currentButtonState));
	memset(previousButtonState, false, sizeof(previousButtonState));

	memset(currentKeyState, false, sizeof(currentKeyState));
	memset(previousKeyState, false, sizeof(previousKeyState));

	scrollDelta = 0.0f;
}

InputHandler::~InputHandler()
{

}

bool InputHandler::update()
{
	for (unsigned short i = 0; i < sf::Keyboard::KeyCount; ++i)
	{
		previousKeyState[i] = currentKeyState[i];
		currentKeyState[i] = sf::Keyboard::isKeyPressed(static_cast<sf::Keyboard::Key>(i));
	}

	for (unsigned short i = 0; i < sf::Mouse::ButtonCount; ++i)
	{
		previousButtonState[i] = currentButtonState[i];
		currentButtonState[i] = sf::Mouse::isButtonPressed(static_cast<sf::Mouse::Button>(i));
	}

	scrollDelta = 0.0f;

	return true;
}
