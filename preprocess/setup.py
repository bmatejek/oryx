from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy as np

extensions = [
    Extension(
        name='components',
        include_dirs=[np.get_include()],
        sources=['components.pyx', 'cpp-components.cpp'],
        extra_compile_args=['-O4', '-std=c++0x'],
        language='c++'
    ),
]

setup(
    name='components',
    ext_modules=cythonize(extensions)
)
