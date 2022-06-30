# Visual Tests

Last run date: 2022-06-30.14:23:31

## Usage 

Visual tests module is developed with the C++17 and requires Eigen3 package:

For Ubuntu it can be installed with 
```
sudo apt install libeigen3-dev
```

To run the visualisation you will need to build the CMake project with the 
VIS=ON and run the executable `./visualtests-out/visualtests` or simply use the makefile command:

```bash 
make build-dev VIS=ON run-vis 
```

<hr>

### RiftingPauline

| Result  | Reference                   |
| ------------- |-----------------------------|
| ![](img/RiftingPauline.png)  | ![](img/RiftingPauline.png) |


### ShearTemplate

| Result  | Reference                   |
| ------------- |-----------------------------|
| ![](img/ShearTemplate.png)  | ![](img/ShearTemplateReference.png) |


### ShearTemplate with shear_style = 1

| Result  | Reference                   |
| ------------- |-----------------------------|
| ![](img/ShearTemplate1.png)  | ![](img/ShearTemplate1Reference.png) |

### ShearHeatingDuretz14

| Result  | Reference                   |
| ------------- |-----------------------------|
| ![](img/ShearHeatingDuretz14.png)  | ![](img/ShearHeatingDuretz14Reference.png) |


### TopoBenchCase1 Result with Analytical solution

<img style="width: 75%" src="img/TopoBenchCase1.png"/>
