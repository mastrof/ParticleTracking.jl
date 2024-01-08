# ParticleTracking.jl

Detection and tracking of blob-like objects in Julia.

The goal of this package is to provide a pure Julia alternative to
libraries available in other languages and frameworks
(e.g. trackpy, TrackMate), aiming at speed and usability
at least for simple tracking problems.

## Features
- Blob detection based on Laplacian of Gaussian filtering
- Subpixel localization
- Min-cost flow tracking with custom cost functions
- Convenient plot recipes and interactive GUIs via Makie.jl

## Acknowledgements
This project has received funding from the European Union's Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie grant agreement No 955910.
