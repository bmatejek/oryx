from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy as np

extensions = [
    Extension(
        name='segment',
        include_dirs=[np.get_include()],
        sources=['segment.pyx', 'cpp-segment.cpp'],
        extra_compile_args=['-O4', '-std=c++0x'],
        language='c++'
    ),
]

setup(
    name='segment',
    ext_modules=cythonize(extensions)
)
