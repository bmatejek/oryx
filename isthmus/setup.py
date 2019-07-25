from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy as np

extensions = [
    Extension(
        name='topological_thinning_downsampled',
        include_dirs=[np.get_include()],
        sources=['topological_thinning_downsampled.pyx', 'cpp-topological-thinning-downsampled.cpp'],
        extra_compile_args=['-O4', '-std=c++0x'],
        language='c++'
    )
]

setup(
    name='topological_thinning_downsampled',
    ext_modules=cythonize(extensions)
)
