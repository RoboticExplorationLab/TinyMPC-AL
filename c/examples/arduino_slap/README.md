# Using `slap` in Arduino Projects

The `slap` library is written to be entirely compatible with
Arduino projects. 

It will (hopefully) be registered as an official Arduino package
at some point, but until then, follow these instructions 
for incorporating it into your Arduino projects.

## Arduino CLI
We recommend using the Arduino CLI, which provides more control
than the Arduino IDE. We have tested this with the Teensy 4.0 
board, which is detailed below.


1. Install Arduino CLI
    ```
   curl -fsSL https://raw.githubusercontent.com/arduino/arduino-cli/master/install.sh | sh
   ```
2. Install Arduino packages
    ```
   arduino-cli core install teensy:avr
   ```
3. Add yourself to the `dialout` group (if needed)
    ```
    sudo usermod -a -G tty <username>
    sudo usermod -a -G dialout <username>
    ```
4. Compile the code
    ```
   cd examples/
   arduino-cli compile -b teensy:avr:teensy40 --libraries ../.. arduino_slap
   ```
   
Can you also add the library to your system Arduino library
path instead of manually specifying the location as we did above.
