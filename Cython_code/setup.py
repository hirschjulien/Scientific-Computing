from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy
import os
import sys

# LAPACK 路径配置
LAPACK_ROOT = "/opt/homebrew/Cellar/lapack/3.12.1"
LAPACK_INCLUDE = os.path.join(LAPACK_ROOT, "include")
LAPACK_LIB = os.path.join(LAPACK_ROOT, "lib")

# Homebrew 路径
BREW_PREFIX = "/opt/homebrew"
gsl_include = f'{BREW_PREFIX}/include'
gsl_lib = f'{BREW_PREFIX}/lib'

extensions = [
    Extension(
        "DP54",
        ["DP54.pyx"],
        include_dirs=[
            numpy.get_include(),
            gsl_include,
            LAPACK_INCLUDE,
        ],
        library_dirs=[
            gsl_lib,
            LAPACK_LIB,
        ],
        libraries=[
            "gsl", 
            "gslcblas",
            "lapack",     # LAPACK库
            "blas",       # BLAS库
            "lapacke"     # LAPACKE库
        ],
        extra_compile_args=[
            "-O2", 
            #"-march=native",           # 添加CPU优化
            #"-ffast-math",            # 启用快速数学运算
            #"-funroll-loops",         # 启用循环展开 
            "-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION",
            f"-I{LAPACK_INCLUDE}"
        ],
        language='c',
        extra_link_args=[
            "-framework", "Accelerate",
            f"-L{LAPACK_LIB}",
            "-llapack",
            "-lblas",
            "-llapacke"
        ],
        runtime_library_dirs=[LAPACK_LIB]  # 运行时库路径
    )
]

setup(
    name="ExplicitEulerSolver",
    ext_modules=cythonize(
        extensions,
        compiler_directives={'language_level': "3", 'profile': False,
                               'linetrace': False,
                               'binding': False},
        annotate=True
    ),
)
