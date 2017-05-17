## Steerable filters
<img style="float: right;" src="http://www.francoisaguet.net/img/dnaHSV.jpg">
Steerable filters for edge and ridge detection, based on [1] and [2]. Implementation in c++, with wrappers for Python and Matlab.

### Python

<!-- Installation:
```
pip install steerablefilter
``` -->

<!-- ```python
import steerablefilter as sf
res, angle = sf.Detector2D()
``` -->

### Matlab

Installation:
1. Add `steerablefilter/matlab` to the Matlab path.
2. Compile the MEX file. [GSL](https://www.gnu.org/software/gsl/) (GNU Scientific Library) must be installed.
```matlab
cd matlab
mex -I/usr/local/include -I../src /usr/local/lib/libgsl.a -output steerableDetector ../src/steerableDetector.cpp ../src/steerableDetector_mex.cpp
```

Documentation:
```matlab
help steerableDetector
```

### References
[1] Jacob and Unser, <i>IEEE Trans. Pattern Anal. Mach. Intell.</i> 26(8), pp. 1007-1019, 2004. [[link](http://ieeexplore.ieee.org/document/1307008/)]</br>
[2] Aguet et al., <i>Proc. IEEE ICIP</i>, pp. II 1158-1161, 2005. [[link](http://ieeexplore.ieee.org/document/1530266/)]
