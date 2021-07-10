#pragma once

#include <SFML/Graphics.hpp>

// translates basic input from keyboard and mouse for more extensible usage
//
class InputHandler
{
public:
	InputHandler();
	~InputHandler();

	// call at start of loop e.g., before poll event
	//
	bool update();

public:
	inline bool left_pressed() const { return (currentButtonState[sf::Mouse::Button::Left] && !previousButtonState[sf::Mouse::Button::Left]); }
	inline bool right_pressed() const { return (currentButtonState[sf::Mouse::Button::Right] && !previousButtonState[sf::Mouse::Button::Right]); }
	inline bool middle_pressed() const { return (currentButtonState[sf::Mouse::Button::Middle] && !previousButtonState[sf::Mouse::Button::Middle]); }

	inline bool left_released() const { return (!currentButtonState[sf::Mouse::Button::Left] && previousButtonState[sf::Mouse::Button::Left]); }
	inline bool right_released() const { return (!currentButtonState[sf::Mouse::Button::Right] && previousButtonState[sf::Mouse::Button::Right]); }
	inline bool middle_released() const { return (!currentButtonState[sf::Mouse::Button::Middle] && previousButtonState[sf::Mouse::Button::Middle]); }

	inline bool left_held() const { return currentButtonState[sf::Mouse::Button::Left]; }
	inline bool right_held() const { return currentButtonState[sf::Mouse::Button::Right]; }
	inline bool middle_held() const { return currentButtonState[sf::Mouse::Button::Middle]; }

	inline bool scroll_up() const { return scrollDelta > 0; }
	inline bool scroll_down() const { return scrollDelta < 0; }

	inline bool key_pressed(sf::Keyboard::Key key) const
	{ 
		if (key < 0 || key >= sf::Keyboard::Key::KeyCount)
			return false;

		return (currentKeyState[key] && !previousKeyState[key]);
	}
	inline bool key_released(sf::Keyboard::Key key) const
	{
		if (key < 0 || key >= sf::Keyboard::Key::KeyCount)
			return false;

		return (!currentKeyState[key] && previousKeyState[key]);
	}
	inline bool key_held(sf::Keyboard::Key key) const
	{
		if (key < 0 || key >= sf::Keyboard::Key::KeyCount)
			return false;

		return currentKeyState[key];
	}

	// call in pollEvent when MouseWheelScrolled event
	//
	void set_scrollDelta(float scrollDelta)
	{
		this->scrollDelta = scrollDelta;
	}

private:
	bool currentButtonState[sf::Mouse::ButtonCount];
	bool previousButtonState[sf::Mouse::ButtonCount];
	
	bool currentKeyState[sf::Keyboard::KeyCount];
	bool previousKeyState[sf::Keyboard::KeyCount];

	float scrollDelta;
};

