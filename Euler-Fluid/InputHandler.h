#pragma once

#include <SFML/Graphics.hpp>

// translates basic input from keyboard and mouse for more extensible usage
//
class InputHandler
{
public:
	InputHandler();
	~InputHandler();

	// call at start of loop before poll event
	//
	bool update();

public:
	inline bool get_left_pressed() const { return (currentButtonState[sf::Mouse::Button::Left] && !previousButtonState[sf::Mouse::Button::Left]); }
	inline bool get_right_pressed() const { return (currentButtonState[sf::Mouse::Button::Right] && !previousButtonState[sf::Mouse::Button::Right]); }
	inline bool get_middle_pressed() const { return (currentButtonState[sf::Mouse::Button::Middle] && !previousButtonState[sf::Mouse::Button::Middle]); }

	inline bool get_left_released() const { return (!currentButtonState[sf::Mouse::Button::Left] && previousButtonState[sf::Mouse::Button::Left]); }
	inline bool get_right_released() const { return (!currentButtonState[sf::Mouse::Button::Right] && previousButtonState[sf::Mouse::Button::Right]); }
	inline bool get_middle_released() const { return (!currentButtonState[sf::Mouse::Button::Middle] && previousButtonState[sf::Mouse::Button::Middle]); }

	inline bool get_left_held() const { return currentButtonState[sf::Mouse::Button::Left]; }
	inline bool get_right_held() const { return currentButtonState[sf::Mouse::Button::Right]; }
	inline bool get_middle_held() const { return currentButtonState[sf::Mouse::Button::Middle]; }

	inline bool get_scroll_up() const { return scrollDelta > 0; }
	inline bool get_scroll_down() const { return scrollDelta < 0; }

	inline bool get_key_pressed(sf::Keyboard::Key key) const
	{ 
		if (key < 0 || key >= sf::Keyboard::Key::KeyCount)
			return false;

		return (currentKeyState[key] && !previousKeyState[key]);
	}

	inline bool get_key_released(sf::Keyboard::Key key) const
	{
		if (key < 0 || key >= sf::Keyboard::Key::KeyCount)
			return false;

		return (!currentKeyState[key] && previousKeyState[key]);
	}

	inline bool get_key_held(sf::Keyboard::Key key) const
	{
		if (key < 0 || key >= sf::Keyboard::Key::KeyCount)
			return false;

		return currentKeyState[key];
	}

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

