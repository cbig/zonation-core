+ [Introduction](https://github.com/cbig/zonation-core#zonation-core)
+ [Authors](https://github.com/cbig/zonation-core#authors-and-contributors)
+ [Compiling from source](https://github.com/cbig/zonation-core#compiling-dynamically-from-source)
+ [Running the tutorial](https://github.com/cbig/zonation-core#testing-zonation-using-the-tutorial-data)
+ [License](https://github.com/cbig/zonation-core#license)
+ [Attribution](https://github.com/cbig/zonation-core#attribution)


## Zonation core

Zonation is a spatial conservation prioritization framework for
large-scale conservation planning. It identifies areas, or landscapes,
important for retaining high habitat quality and connectivity for
multiple biodiversity features (eg. species), providing a quantitative
method for enhancing species' long term persistence.

Essentially, this software is a decision support tool for all
non-commercial parties working on conservation issues. As Zonation
operates on large grids, it provides a direct link between GIS,
statistical distribution modeling and spatial conservation
prioritization.

This repository contains the computational core of Zonation software.
It does not contain the graphical user interface (GUI) components.

You can dowload the latest stable release installers (including the GUI) from [the official website](https://www.helsinki.fi/en/researchgroups/metapopulation-research-centre/software#section-14300).

## Authors and contributors

Authors:

+ Atte Moilanen <<atte.moilanen@helsink.fi>>
+ Federico Montesino Pouzols
+ Jarno Leppänen

Contributors:

+ Joona Lehtomäki <<joona.lehtomaki@helsink.fi>>

## Compiling (dynamically) from source

These instruction describe the dynamic compilation of Zonation 4.0.0 on Ubuntu 14.04 (Trusty Tahr). Instructions have 
not been tested on other distributions/versions. For a set of bash scripts automating the steps below, please see 
[this repository](https://github.com/cbig/zig4-compilation-scripts). 

#### 1. Install Zonation dependencies

```
sudo apt-get update

sudo apt-get install cmake build-essential unzip libqt4-dev libfftw3-dev libqwt-dev libboost-all-dev libgdal-dev
``` 

#### 2. Get Zonation sources

Create a suitable download directory, fetch Zonation source code (version 4.0.0) from CBIG server, and extract the sources to the directory created:

```
mkdir zonation
wget https://github.com/cbig/zonation-core/archive/master.zip -P zonation
unzip zonation/master.zip -d zonation
```

or if you're using git, clone this repo 

```
mkdir zonation
git clone https://github.com/cbig/zonation-core.git zonation/zonation-core-master
```

#### 3. Build Zonation

To build Zonation library (`zig4lib`) and Zonation CLI utility (`zig4`), do the following:

```
mkdir zonation/build
cd zonation/build
cmake ../zonation-core-master
make
```

If you have several cores available for compilation, you can pass switch `-jX` to `make` where `X` is the number of 
designated cores (e.g. `make -j4`).

#### 4. Make Zonation available system wide (optional)

In order to call `zig4` anywhere on the system (instead of just the build-location), create a symbolic link:

```
sudo ln -s FULL_PATH/zonation/build/zig4/zig4 /usr/local/bin/zig4
```

Replace `FULL_PATH` with the full path to the directory containing directory `zonation` created in step 2.

## Runnning Zonation using the tutorial data

This stage is completely optional.

You can test that Zonation works by running some of the runs used in the 
[Hunter Valley tutorial](https://github.com/cbig/zonation-tutorial). This is not a thorough test, but it will at least give you an overview on whether Zonationis working as intended. Tutorial data and setup files are fetched using 
[git](http://git-scm.com/) and run using `zrunner` command line utility found in package [ztools](https://github.com/cbig/ztools).

#### 5. Install ztools dependencies

First, install git:

```
sudo apt-get -y install git
```

Then, install Python packages needed by zrunner:

```
sudo apt-get -y install python-yaml python-pip 
```

#### 6. Install ztools

`ztools` is installed directly from GitHub using [pip](http://www.pip-installer.org/en/latest/).

```
sudo pip install https://github.com/cbig/ztools/archive/master.zip
```

#### 7. Clone Zonation tutorial using git

With git installed, clone the tutorial repository with the following command:

```
git clone https://github.com/cbig/zonation-tutorial.git
``` 

#### 8. Run the tutorial runs

You can run [5 basic tutorial variants](https://github.com/cbig/zonation-tutorial/tree/master/basic) defined in the configuration file `basic/tests/ztests_basic.yaml ` by using zrunner:

```
zrunner -l tutorial_runs.yaml
```

zrunner will produce an output file `results_XXX.yaml` in the same folder. `XXX` will correspond to information about your system. If everything went fine, you should see no critical errors on the screen and the yaml-file should report execution times for successful runs.

----

## License

Zonation computational core (zig4) is distributed under the 
GNU General Public License (GPL) version 3
(http://www.gnu.org/licenses/). 

Zonation computational core
© 2011-2014 Conservation Biology Informatics Group
© 2004-2011 Atte Moilanen

Zonation computational core is free software: you can redistribute it
and/or modify it under the terms of the GNU General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Do not use this software if you disagree with the disclaimer or
conditions of use. Even though the Zonation software has been done
with the best of intentions, it is quite beyond one small research
group to ensure its correct operation under all operating systems and
environments. Anticipating all potential combinations of erroneous
input has not been possible. Therefore, use the software with care and
make an effort to understand how the inputs connect to outputs.

## Attribution

Zonation uses the following libraries: 
- boost (http://www.boost.org)
- FFTW (http://www.fftw.org). We also used
- GDAL (http://www.gdal.org)
- Qt (http://qt-project.org)
- Qwt (http://qwt.sourceforge.net)

Various versions of GCC, the GNU compiler collection (http://gcc.gnu.org) 
were used as well.

These free open source software projects are greatly acknowledged!
