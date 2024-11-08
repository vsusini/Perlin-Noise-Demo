# Perlin Noise Terrain Demo

[![Live Demo](https://img.shields.io/badge/demo-online-green.svg)](https://vsusini.github.io/Perlin-Noise-Demo)

An interactive 3D terrain generator using Perlin and Simplex noise, built with Three.js. This project visualizes noise algorithms through a dynamic point-based terrain system with customizable parameters.

## Overview

This project was developed as part of a fourth-year Computer Graphics course, exploring the implementation and applications of Perlin noise. The visualization uses Three.js to create a terrain-like surface where dot heights are determined by noise functions, with color variations based on elevation.

## Features

- Real-time terrain generation using both Perlin and Simplex noise
- Interactive controls for:
  - Terrain resolution (rows/columns)
  - Wave height and scale
  - Dot size
  - Color customization for different elevation levels
- Multiple camera angles:
  - Sky view
  - Wave-height view
  - First-person perspective
- Dynamic seed generation
- Responsive design

## Technologies Used

- Three.js
- dat.GUI
- HTML5/CSS3
- JavaScript

## Local Development

1. Clone the repository:
   ```bash
   git clone https://github.com/vsusini/Perlin-Noise-Demo.git
   ```

2. Open `index.html` in your preferred browser

3. Start experimenting with the controls in the GUI panel

## Controls

The demo includes several controllable parameters:

### Wave Settings
- Row/Column Number: Adjust terrain resolution
- Wave Height: Control the amplitude of the terrain
- Perlin Scale: Adjust the noise frequency
- Dot Size: Change the size of individual points
- Generate Seed: Create new random terrain patterns
- Noise Algorithm: Toggle between Perlin and Simplex noise

### Color Settings
- Customize colors for different elevation levels:
  - Top Color
  - Top-Middle Color
  - Middle Color
  - Floor Color

### Camera Options
- Sky Camera: Bird's eye view
- Wave-Height Camera: Mid-level perspective
- First Person: Immersive ground-level view

## Acknowledgments

- Based on Stefan Gustavson's Perlin noise implementation
- Optimizations by Peter Eastman
- JavaScript conversion by Joseph Gentle

## Author

Vincent Susini - [Portfolio](https://www.vsusini.com/)

